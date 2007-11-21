/****************** Flame2DState.h **************************************
  This class defines the state variables and constructors for the 
  2D Axisymmertric Navier-Stokes with Multiple Species solution.

  NOTES:  
         - be careful with set_species_data() and set_initial_values(),
           due to dynamic memory 

         - if #define STATIC_NUMBER_OF_SPECIES is set dynamic memory
           is not used so code is faster, however code is not as flexible
           as only up to STATIC_NUMBER_OF_SPECIES can be used without recompling.
  
         - Avoid cState.T() as much as possible as it is $$$ to calculate, pass
           the temperature from pState if possible as it is cheap (p/rhoR)             
              
         - The Static "specdata" has to be deleted outside the class
           in order to avoid seg faults. In Flame2D it is done by calling
           the deallocate function from the Input destructor as it is kept
           until the end of the program.

  TODO:  - Unphysical_properties & negativespec_check functions need to be  
           worked through

                   
***********************************************************************/
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
#include "../Physics/Species.h"
#include "../Physics/NASAData/NASARP1311data.h"
#include "../Reactions/Reactions.h"


//Temperature convergence tolerance in
//Flame2D_cState::T(void)
// these should be moved to CFD.h or Math.h
#define CONV_TOLERANCE  1e-8   //Tolerance used for temperature convergence
#define SPEC_TOLERANCE  1e-8   //Used in negative_speccheck for species round off (was MICRO)
#define TLOWERBOUNDS   200.0   //Uncoment this fix a lower bounds in T calculation

//number of fixed variables in the Flame2D class
#define NUM_FLAME2D_VAR_SANS_SPECIES 4  //rho, v(2), p

// If you define this variable, the number of species will be
// predetermined for faster calculations.., however it is not as general 
#define STATIC_NUMBER_OF_SPECIES 36 //2 AIR, 6 2STEP_CH4


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
class Flame2D_pState {

private:
  //@{ @name Temporary storage:
  static double *vis; // storage for component viscosity
  static double *e_k; // storage for component energy
  //@}
  
  //all public ....
protected:
public:
  //@{ @name Primitive variables and associated constants:
  double         rho;   //!< Density.
  Vector2D         v;   //!< Flow velocity (2D)  
  double           p;   //!< Pressure.
#ifdef STATIC_NUMBER_OF_SPECIES
  Species        spec[STATIC_NUMBER_OF_SPECIES];
#else 
  Species      *spec;   //!< Species class c[ns]
#endif
  double           k;   //!< Turbulent kinetic energy.
  double       omega;   //!< Turbulent specific dissipation rate.
  Tensor2D                    tau; //!< Shear Stress Tensor
  Vector2D                  qflux; //!< Heat Flux Vector  
  
  
  //! Static Variaables 
  static int                   ns; //!< number of species
  static NASARP1311data *specdata; //!< Global Species Data
  static double          *Schmidt; //!< Schmidt Number for each species
  static Reaction_set       React; //!< Global Reaction Data
  static double    low_temp_range; //!< Low temp data range
  static double   high_temp_range; //!< High temp data range
  static int       NUM_VAR_FLAME2D; //!< Number of Chem2d variables (6+ns)
  static double              Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static double         gravity_z; //!< m/s^2 acceleration due to gravity  
  //@}

  //default constructors of many flavours, hopefully one is right 4U   
  //v.zero(); tau.zero(); qflux.zero(); lambda.zero(); theta.zero();
  //  ^ Unecessary as Vector2D and Tensor2D defaults are ZERO
  //@{ @name Creation, copy, and assignment constructors.
 
  Flame2D_pState(): rho(DENSITY_STDATM), p(PRESSURE_STDATM)
  { specnull();  set_initial_values(); }
  
  Flame2D_pState(const double &value): rho(value), v(value), p(value)
  { specnull(); set_initial_values(value); }
  
  Flame2D_pState(const double &d, const Vector2D &V, 
		const double &pre):
    rho(d), v(V), p(pre)
  { specnull();  set_initial_values(); }
  
  Flame2D_pState(const double &d, const double &vx, 
		const double &vy, const double &pre):
    rho(d), v(vx,vy), p(pre)
  {specnull();  set_initial_values(); }

  Flame2D_pState(const double &d, const Vector2D &V, 
		const double &pre, const double &value):
    rho(d), v(V), p(pre)
  { specnull();  set_initial_values(value); }

  Flame2D_pState(const double &d, const double &vx, 
		const double &vy, const double &pre, 
		const double &value):
    rho(d), v(vx,vy), p(pre)
  { specnull();  set_initial_values(value); }
  
  Flame2D_pState(const double &d, const double &vx, 
		const double &vy, const double &pre, 
		const Species *mfrac): 
    rho(d), v(vx,vy), p(pre)
  { specnull();  set_initial_values(mfrac); }
  
  Flame2D_pState(const double &d, const Vector2D &V, 
		const double &pre, const Species *mfrac): 
    rho(d), v(V), p(pre)
  { specnull(); set_initial_values(mfrac); }
  
  //this is needed for the operator overload returns!!!!
  Flame2D_pState(const Flame2D_pState &W): rho(W.rho), v(W.v), p(W.p),
 					 tau(W.tau), qflux(W.qflux)
  { specnull(); set_initial_values(W.spec); }      
  //@}

   //read in ns species data, call only once as its static
   void set_species_data(const int &,const string *,const char *,
 			 const double&,const double *, const int &);

   //set initial data values predominately used internally !!!
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *cfrac);
   void set_initial_values(const Species *mfrac);

   //Copy construtor, cheaper than = operator
   void Copy(const Flame2D_pState &W);

   // return the number of variables - number of species
   int NumVarSansSpecies() const { return NUM_FLAME2D_VAR_SANS_SPECIES; }

   /*************** VACUUM OPERATOR *********************/
  void Vacuum(){ rho=ZERO; v.zero(); p=ZERO;
     for(int i=0; i<ns; i++)  spec[i].Vacuum();
     tau.zero(); qflux.zero();
   }

   void zero_non_sol(){
     for(int i=0; i<ns; i++){
       spec[i].gradc.zero();
       spec[i].diffusion_coef=ZERO;
     }
     tau.zero(); qflux.zero();
   }  
 
  //! set acceleration due to gravity
  void set_gravity(const double &g);
  //@}

   /******** Set Data Temperature Ranges ***************/
   void Temp_low_range();     
   void Temp_high_range(); 

   /***************** Mixing Rules ********************
    The following constructors return "total" physical
    parameters based on mixture rules for each.
   ****************************************************/
   double Mass(void) const;   //mixture molecular mass
   double Rtot(void);         //mixture gas constant
   double Rtot(void) const;   
   double Cp(void) const;     //mixture heat capacity (Pressure constant)
   double Cp(const double& TEMP) const;
   double Cv(void) const;     //mixture heat capacity (Volume constant)
   double g(void) const;      //mixture heat capacity ratio
   double e(void) const;      //mixture absolute (sensible+chemical) internal energy
   double eref(void) const;
   double es(void) const;     //mixture sensible internal energy  
   double h(void) const;      //mixture specific enthalpy
   double h(const double &T) const;
   double href(void) const;
   double hs(void) const;
   double hs(const double &T) const;
   double E(void) const;      //mixture total internal energy  
   double H(void) const;      //mixture total enthalpy
   double Hs(void) const;     //mixture total enthalpy
   double mu(void) const;     //mixture viscosity
   double kappa(void) const;  //mixture thermal conductivity
   double hprime(void) const;  
   double hprime(double &Temp) const;
   //double Diffusion_coef() const;  //mixture diffusion coefficient 

   /***************************************************/
   Vector2D rhov(void) const;      //Momentum
   double T(void) const;           //Temperature
   double T(double &h_s) const;    //Determine temperature knowing sensible enthalpy
   double gamma_guess(void) const; //gamma for iterative schemes
   double a(void);                 //speed of sound
   double a(void) const;
//    bool negative_speccheck(void) ; //-ve mass frac check and sets small -ve c's to ZERO  
//    bool Unphysical_Properties_Check(Flame2D_cState &U, const int Flow_Type, const int n) ;   
   double SpecCon(int i) const;     //Species i concentration (rho*c/mol_mass)
   double MoleFrac(int i) const;    //Species i mole fraction (c*mol_mass_mix/mol_mass)
   double Gibbs(int species) const; //Gibbs Free Energy (H-TS) for species

   /********** radiating species related parameters ***********/
    void MoleFracOfRadSpec( double &xCO,  double &xH2O, 
			    double &xCO2, double &xO2 ) const;

   /************* Dimensionless Parameters ********************/
   double Schmidt_No(const int &) const;
   double Prandtl() const;
   double Lewis(const int &) const;

   /************** Temperature Derivatives *******************/
   double diedip() const;
   double diedirho() const;
   double dmudT(void) const;
   double dkappadT(void) const;

   /************ Strain rate tensor, laminar stress and Reynolds stress tensors ***********/
   Tensor2D Strain_Rate(const Flame2D_pState &dWdx,
			const Flame2D_pState &dWdy,
			const int Axisymmetric,
			const Vector2D X);

   void Laminar_Stress(const Flame2D_pState &dWdx,
		       const Flame2D_pState &dWdy,
		       const int Axisymmetric,
		       const Vector2D X);

   /************ Heat Flux vector thermal Diffusion ***********/
   Vector2D thermal_diffusion(void) const;

   /*************** Conserved solution state. ****************/
   Flame2D_cState U(void) const;
   Flame2D_cState U(const Flame2D_pState &W) const;
   friend Flame2D_cState U(const Flame2D_pState &W);

   /**************** Fluxes ***********************************/
   /* Inviscid Solution flux (x-direction). */
   Flame2D_cState Fx(void) const;

   /*********** Inviscid Flux Jacobian X ***********************/
   friend void dFIdU(DenseMatrix &dFdU, const Flame2D_pState &W);
   friend void dFIdU_FD(DenseMatrix &dFdU, const Flame2D_pState &W);
   friend void dFIdW(DenseMatrix &dFdW, const Flame2D_pState &W);
   friend void dFIdW_FD(DenseMatrix &dFdW, const Flame2D_pState &W);

   /** Conserved/Primitive  & Primitive/Conservered Jacobians **/
   void dWdU(DenseMatrix &dWdQ)const;
   void dWdU_FD(DenseMatrix &dWdQ);

   void dUdW(DenseMatrix &dQdW);
   void dUdW_FD(DenseMatrix &dUdW);

   /************* Eigenvalues *********************************/
   /* Eigenvalue(s) (x-direction). */ 
   Flame2D_pState lambda_x(void) const;

   /************** Precondtioned Eigenvalues ******************/
   Flame2D_pState lambda_preconditioned_x(const double &MR2) const;

   /**************** Eigenvectors *****************************/  
   Flame2D_cState rc_x(const int &index) const; // Conserved right (x-direction)
   Flame2D_pState lp_x(const int &index) const; // Primitive left (x-direction)

   /************** Preconditioned Eigenvectors ****************/
   Flame2D_cState rc_x_precon(const int &index,const double &MR2) const;  // Conserved right (x-direction)
   Flame2D_pState lp_x_precon(const int &index,const double &MR2) const;  // Primitive left  (x-direction)

   /*************** Preconditioner ****************************/
   double u_plus_aprecon(const double &u,const int &flow_type_flag,const double &deltax) const;
   void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
   double Mr2(const int &flow_type_flag, const double &deltax) const;
   void Low_Mach_Number_Preconditioner(DenseMatrix &P,const int &flow_type_flag, const double &deltax) const; 
   void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,const int &flow_type_flag, const double &deltax) const; 

   /********** Axisymmetric Source Terms **********************/
   Flame2D_cState Sa_inviscid(const Vector2D &X,
                             const int Axisymmetric) const;

   Flame2D_cState Sa_viscous(const Flame2D_pState &dWdx, 
                            const Flame2D_pState &dWdy, 
                            const Vector2D &X,
                            const int Axisymmetric);

   //Due to axisymmetric coordinate system
  void dSa_idU(DenseMatrix &dSa_idU, const Vector2D &X, const int Axisymmetric) const;  //Inviscid source Jacobian
  void dSa_idU_FD(DenseMatrix &dSa_idU, const Vector2D &X,const int Axisymmetric) const;
  void dSa_vdW(DenseMatrix &dSa_VdW, const Flame2D_pState &dWdx,
 	       const Flame2D_pState &dWdy,const Vector2D &X, 
 	       const int Axisymmetric, 
 	       const double d_dWdx_dW, const double d_dWdy_dW) const;  //Viscous source Jacobian

   /****** Source terms associated with finite-rate chemistry ******/
   Flame2D_cState Sw(int &REACT_SET_FLAG ) const;
   void dSwdU(DenseMatrix &dSwdU,const int &solver_type) const; //Jacobian
   void dSwdU_FD(DenseMatrix &dSwdU) const;
   double dSwdU_max_diagonal(const int &Preconditioned,					 
 			    const double &delta_n,
 			    const int &solver_type) const;

   /****** Source terms associated with gravitational forces ******/
   Flame2D_cState Sg(void) const;
   void dSgdU(DenseMatrix &dSgdU) const; //Jacobian

  /****** Source terms associated with dual time stepping *******/
  Flame2D_cState S_dual_time_stepping(const Flame2D_cState &U,
                                     const Flame2D_cState &Ut,
                                     const Flame2D_cState &Uold,
                                     const double &dTime,
                                     const int &first_step) const;

   /**************** Operators Overloading ********************/
   /* Index operator */
   double &operator[](int index);
   const double &operator[](int index) const;

   /* Binary arithmetic operators. */
   Flame2D_pState operator +(const Flame2D_pState &W) const;
   Flame2D_pState operator -(const Flame2D_pState &W) const;
   Flame2D_pState operator *(const double &a) const;
   friend Flame2D_pState operator *(const double &a, const Flame2D_pState &W);
   Flame2D_pState operator /(const double &a) const;
   double operator *(const Flame2D_pState &W) const;
   Flame2D_pState operator ^(const Flame2D_pState &W) const;

   /* Assignment Operator. */ 
   Flame2D_pState& operator =(const Flame2D_pState &W); 

   /* Shortcut arithmetic operators. */
   Flame2D_pState& operator +=(const Flame2D_pState &W);
   Flame2D_pState& operator -=(const Flame2D_pState &W);
 
   /* Unary arithmetic operators. */
   //Flame2D_pState operator +(const Flame2D_pState &W);
   friend Flame2D_pState operator -(const Flame2D_pState &W);

   /* Relational operators. */
   friend int operator ==(const Flame2D_pState &W1, const Flame2D_pState &W2);
   friend int operator !=(const Flame2D_pState &W1, const Flame2D_pState &W2);

   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, const Flame2D_pState &W);
   friend istream& operator >> (istream &in_file,  Flame2D_pState &W);

   /**************** Destructors ******************************/
   //for static specdata
   void Deallocate_static(void){ 
     if(specdata != NULL) { delete[] specdata; specdata = NULL; }
     if(Schmidt != NULL) { delete[] Schmidt; Schmidt = NULL; }
     if(vis != NULL) { delete[] vis; vis = NULL; }
     if(e_k != NULL) { delete[] e_k; e_k = NULL; }
   }

#ifdef STATIC_NUMBER_OF_SPECIES    
   void specnull() {}
   void spec_memory() {}
   ~Flame2D_pState(){}
#else
   void specnull() {spec=NULL;}
   void spec_memory() { Deallocate(); spec = new Species[ns];}
   void Deallocate(void){ if(spec != NULL){ delete[] spec;}  specnull();  }
   ~Flame2D_pState(){ Deallocate(); }
#endif
 
   /**************** Species Resize ***************************/
  void resize_species();
  void resize_species(const int n_old, const string* species_old);
};


 /**************************************************************************
 ********************* FLAME2D_CSTATE CLASS DEFINTION ***********************
     The conserved variables edition of FLAME2D_PSTATE... 
 ***************************************************************************
 ***************************************************************************/
 class Flame2D_cState {
   private: 
   //@{ @name Temporary storage:
   static double *vis; // storage for component viscosity
   static double *e_k; // storage for component energy
   //@}

   //all public .... yes I know ....
   protected:
   public:
   /****** Conservative Vector "U" **************************/
  double                      rho; //!< Density.
  Vector2D                   rhov; //!< Momentum (2D)
  double                        E; //!< Total Energy (rho *(e + HALF*v^2))
#ifdef STATIC_NUMBER_OF_SPECIES
   Species   rhospec[STATIC_NUMBER_OF_SPECIES];
#else 
  Species                *rhospec; //!< Species class using (rho*c[ns])
#endif
   Tensor2D                    tau; //!< Shear Stress Tensor
   Vector2D                  qflux; //!< Heat Flux Vector  

   static int                   ns; //!< number of species
   static NASARP1311data *specdata; //!< Global species data 
   static double          *Schmidt; //!< Schmidt Number for each species
   static double    low_temp_range; //!< Low temp data range
   static double   high_temp_range; //!< High temp data range
   static int       NUM_VAR_FLAME2D; //!< Number of Chem2d variables (4+ns)
   static double              Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
   static double         gravity_z; //!< m/s^2 acceleration due to gravity  
   //@}
 
   /*****************************************************/
   //default constructors of many flavours, hopefully one is right 4U
   Flame2D_cState(): rho(DENSITY_STDATM), E(PRESSURE_STDATM/(rho*((double)0.4)))
   {rhospecnull(); set_initial_values(); }   

   Flame2D_cState(const double &value): rho(value), rhov(value), E(value)
   { rhospecnull();  set_initial_values(value); }   

   Flame2D_cState(const double &d, const Vector2D &V, 
		 const double &En):
     rho(d), rhov(V), E(En)
   { rhospecnull(); set_initial_values(); }

   Flame2D_cState(const double &d, const Vector2D &V, 
		 const double &En, const double &value):
     rho(d), rhov(V), E(En)
   { rhospecnull(); set_initial_values(value); }

   Flame2D_cState(const double &d, const double &vx, 
		 const double &vy, const double &En):
     rho(d), rhov(vx,vy), E(En)
   { rhospecnull(); set_initial_values(); }

   Flame2D_cState(const double &d, const double &vx, 
		 const double &vy, const double &En, 
		 const double &value): 
     rho(d), rhov(vx,vy), E(En)
   { rhospecnull();  set_initial_values(value); }

   Flame2D_cState(const double &d, const double &vx, 
		 const double &vy, const double &En,   
		 const Species *rhomfrac): 
     rho(d), rhov(vx,vy), E(En)
   { rhospecnull();   set_initial_values(rhomfrac); }

   Flame2D_cState(const double d, const Vector2D &V, 
		 const double &En, const Species *rhomfrac): 
     rho(d), rhov(V), E(En)
   { rhospecnull();   set_initial_values(rhomfrac); }
		 
   // WARNING - automatic type conversion
   Flame2D_cState(const Flame2D_pState &W) :  rho(W.rho), rhov(W.rhov()),
					    E(W.E()),
					    tau(W.tau), qflux(W.qflux)
   {   
     for(int i=0; i<W.ns; i++){
       rhospec[i].c = W.rho*W.spec[i].c;
       rhospec[i].gradc = W.rho*W.spec[i].gradc;
       rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
     }  
   }



   //this is needed for the operator overload returns!!!!
   Flame2D_cState(const Flame2D_cState &U): rho(U.rho), rhov(U.rhov), E(U.E),
					  tau(U.tau), qflux(U.qflux)
   { rhospecnull(); set_initial_values(U.rhospec); }

   //read in ns species data, call only once as its static
   void set_species_data(const int &,const string *,const char *,
			 const double&,const double *, const int &);

   //set initial data values predominately used internally 
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *rhomfrac);
   void set_initial_values(const Species *rhomfrac);

   //Copy construtor, cheaper than = operator
   void Copy(const Flame2D_cState &U);

   /***************** VACUUM ************************/
   void Vacuum(){ rho=ZERO; rhov.zero(); E=ZERO;
     for(int i=0; i<ns; i++) rhospec[i].Vacuum();
     tau.zero();  qflux.zero();
   }  

   void zero_non_sol(){
     for(int i=0; i<ns; i++){
       rhospec[i].gradc.zero();
       rhospec[i].diffusion_coef=ZERO;
     }
     tau.zero(); qflux.zero();
   }  

   //! set acceleration due to gravity
   void set_gravity(const double &g);
   //@}

   /******** Set Data Temperature Ranges ***************/
   void Temp_low_range();     
   void Temp_high_range(); 

   /***************** Mixing Rules ********************
    The following constructors return "total" physical
    parameters based on mixture rules for each.
   ****************************************************/
   //double Mass(void);       //mixture molecular mass
   double Rtot(void) const; 
 //  double Cp(void) const;   //mixture heat capacity (Pressure constant)
 //  double Cv(void) const;   //mixture heat capacity (Volume constant)
 //  double g(void) const;    //specific heat ratio
   double gamma_guess(void) const;   //mixture specifc heat ratio
   double e(void) const;             //mixture specific internal energy
   double es(void) const;            //mixture specific sensible internal energy
   double h(const double &T) const;  //mixture specific enthalpy
   double hs(const double &T) const; //mixture specific enthalpy
   double hprime(const double &T) const; 
   double heatofform(void) const;      //mixture heat of formations
   //double E(void) const;             //mixture total internal energy
   //double H(void);                   //mixture total enthalpy
   double mu(void) const;              //mixture viscosity
   double kappa(void) const;           //mixture thermal conductivity

   /***********************************************************/
   Vector2D v(void) const;    //velocity  
   double p(void) const;      //pressure
   double T(void) const;      //temperature
   double a(void) const;      //speed of sound
   bool negative_speccheck(const int &step) ; //-ve mass frac check and sets small -ve c's to ZERO
   bool Unphysical_Properties();
   bool Unphysical_Properties_Check(const int n);   
   double sum_species(void) const;

   //@{ @name Functions required for multigrid.
   //! Copy variables solved by multigrid only.
   void Copy_Multigrid_State_Variables(const Flame2D_cState &Ufine);
  
   //! Zero variables not-solved by multigrid.
   void Zero_Non_Multigrid_State_Variables(void);
   //@}


   /************** Temperature Derivatives *******************/
   double dmudT(void) const;

   // STRAIN RATE, LAMINAR STRESS, REYNOLD STRESS 6 FUNCTIONS\
   // LIKE PRIMITIVE 

   /************ Heat Flux vector thermal Diffusion ***********/
   Vector2D thermal_diffusion(const double &Temp) const;

   /*********** Primitive solution state ***********************/
   Flame2D_pState W(void) const;
   Flame2D_pState W(const Flame2D_cState &U) const;
   friend Flame2D_pState W(const Flame2D_cState &U);

   /**************** Fluxes ***********************************/
   //Viscous Flux (laminar+turbulent)
   Flame2D_cState Viscous_Flux_x(const Flame2D_pState &dWdx) const;
   Flame2D_cState Viscous_Flux_y(const Flame2D_pState &dWdy) const;

   /*************** Preconditioner ****************************/
   void Low_Mach_Number_Preconditioner(DenseMatrix &P,const int &flow_type_flag, const double &deltax) const; 
   void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,const int &flow_type_flag, const double &deltax) const; 

   /***************** Index operators *************************/
   double &operator[](int index);
   const double &operator[](int index) const;

   /**************** Operators Overloading ********************/
   /* Binary arithmetic operators. */
   Flame2D_cState operator +(const Flame2D_cState &U) const;
   Flame2D_cState operator -(const Flame2D_cState &U) const;
   Flame2D_cState operator *(const double &a) const;
   friend Flame2D_cState operator *(const double &a, const Flame2D_cState &U);
   Flame2D_cState operator /(const double &a) const;

   double operator *(const Flame2D_cState &U) const;
   Flame2D_cState operator ^(const Flame2D_cState &U) const;

   /* Assignment Operator. */ 
   Flame2D_cState& operator =(const Flame2D_cState &U); 

   /* Shortcut arithmetic operators. */
   Flame2D_cState& operator +=(const Flame2D_cState &U);
   Flame2D_cState& operator -=(const Flame2D_cState &U);
   Flame2D_cState& operator *=(const double &a);
   Flame2D_cState& operator /=(const double &a);

   /* Unary arithmetic operators. */
   //Flame2D_cState operator +(const Flame2D_cState &U);
   friend Flame2D_cState operator -(const Flame2D_cState &U);

   /* Relational operators. */
   friend int operator ==(const Flame2D_cState &U1, const Flame2D_cState &U2);
   friend int operator !=(const Flame2D_cState &U1, const Flame2D_cState &U2);

   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, const Flame2D_cState &U);
   friend istream& operator >> (istream &in_file,  Flame2D_cState &U);

   /**************** Destructors ******************************/
   void Deallocate_static(void){ 
     if(specdata != NULL) { delete[] specdata; specdata = NULL; }
     if(Schmidt != NULL)  { delete[] Schmidt;   Schmidt = NULL; }
     if(vis != NULL) { delete[] vis; vis = NULL; }
     if(e_k != NULL) { delete[] e_k; e_k = NULL; }
   }

#ifdef STATIC_NUMBER_OF_SPECIES  
   void rhospecnull() {}          
   void rhospec_memory() {}
   ~Flame2D_cState() {}
#else
   void rhospecnull() {rhospec=NULL;}
   void rhospec_memory() { Deallocate(); rhospec = new Species[ns];}
   void Deallocate(void){ if(rhospec != NULL) { delete[] rhospec;} rhospecnull(); }
   ~Flame2D_cState() { Deallocate(); }
#endif

   /**************** Species Resize ***************************/
   void resize_species();
   void resize_species(const int n_old, const string* species_old);

 };

 /**************************************************************************
 ********************* FLAME2D_PSTATE CONSTRUCTORS **************************
 ***************************************************************************/

 /**************************************************************************
    Set up mass fractions memory and initial values, ie. Species Class  
 **************************************************************************/
 inline void Flame2D_pState::set_initial_values(){
   spec_memory();
   for(int i=0; i<ns; i++) spec[i].c = ONE/ns; 
 }

 inline void  Flame2D_pState::set_initial_values(const double &value){
   spec_memory();
   for(int i=0; i<ns; i++) spec[i].c = value; 
 }

 //user specified
 inline void Flame2D_pState::set_initial_values(double *cfrac){
   spec_memory();
   for(int i=0; i<ns; i++) spec[i].c = cfrac[i];
 }

 //another set using species class
 inline void Flame2D_pState::set_initial_values(const Species *mfrac){
   spec_memory();
   for(int i=0; i<ns; i++) spec[i] = mfrac[i];
 }

/**************************************************************************
  Resize the species array, and just set the new array to zero
**************************************************************************/
inline void Flame2D_pState::resize_species() {

  // the new number of species
  int n_new( React.num_species );

  // deallocate first
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES < n_new ) { 
    cerr<<"\n WARNING USING STATIC FLAME2D BUILT WITH "<<STATIC_NUMBER_OF_SPECIES 
	<<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<<n_new<<endl; 
    exit(1); 
  }
#else 
  Deallocate();
#endif  
  
  // set the new value
  ns = n_new;
  
  // reallocate species array and set to zero
  spec_memory();
  set_initial_values(ZERO);
  
}

/**************************************************************************
  Resize the species array, keeping the values
**************************************************************************/
inline void Flame2D_pState::resize_species( const int n_old,           // the old numbre of species
					   const string* species_old) // the old species names array
{
  // temporary storage
  double* c_old = new double[n_old];
  int index;

  // the new number of species
  int n_new( React.num_species );

  // store the old values
  for (int k=0; k<n_old; k++) c_old[k] = spec[k].c;

  // deallocate
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES < n_new ) { 
    cerr<<"\n WARNING USING STATIC FLAME2D BUILT WITH "<<STATIC_NUMBER_OF_SPECIES 
	<<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<<n_new<<endl; 
    exit(1); 
  }
#else 
  Deallocate();
#endif  
  
  // set the new value
  ns = n_new;
  
  // reallocate species array and set to zero
  spec_memory();
  set_initial_values(ZERO);
    
  // match up the species and set their respective values
  for (int k=0; k<n_old; k++) {
    index = React.SpeciesIndex(species_old[k]);
    if (index>=0) spec[index].c = c_old[k];
    else {
      cerr << "\n Flame2D_pState::resize_species() - Species " 
	   << species_old[k]
	   << " not found in current reaction mechanism.\n";
      exit(-1);
    }
  }// endfor - species
  
  // delete temporary memory
  delete[] c_old;

}


 /**********************************************************************
  * Flame2D_pState::set_gravity -- Set the acceleration due to gravity  *
  *                               in m/s^2.  It acts downwards in the  *
  *                               z-dir (g <= 0)                       *
  **********************************************************************/
inline void Flame2D_pState::set_gravity(const double &g) { // [m/s^2]

  // if gravity is acting upwards (wrong way)
  if (g>0) {
    cerr<<"\n Flame2D_pState::set_gravity() - Gravity acting upwards!!!! \n";
    exit(1);
    
  // gravity acting downwards (-ve), OK
  } else {
    gravity_z = g;
  }
}


 /*****************  Momentum *******************************/
 inline Vector2D Flame2D_pState::rhov(void) const{
   return rho*v; 
 } 

 /********************** Prandtl ****************************/
 inline double Flame2D_pState::Prandtl(void) const{
   //Pr = Cp*mu/k
   return Cp()*mu()/kappa();
 }

 /********************** Schmidt ****************************/
 inline double Flame2D_pState::Schmidt_No(const int &i) const{

   if(spec[i].diffusion_coef > ZERO){
     return mu()/(rho*spec[i].diffusion_coef);
   } else {
     return Schmidt[i];
   }
 }

 /********************** Lewis *****************************/
 inline double Flame2D_pState::Lewis(const int &i) const{

   if(spec[i].diffusion_coef > ZERO){
     return kappa()/(rho*Cp()*spec[i].diffusion_coef);
   }
   return ZERO;
 }

 /******* Mixture Diffusion Coefficient ********************/
 // inline double Flame2D_pState::Diffusion_coef(void) const{
 //   double sum=ZERO;
 //   for(int i=0; i<ns; i++){
 //     sum += spec[i].c * spec[i].diffusion_coef;
 //   }
 //   return sum;
 // }

 /************* Temperature ********************************/
 inline double Flame2D_pState::T(void) const{
   return p/(rho*Rtot());
 }

 /************* Strain rate tensor ***********************************/
 inline Tensor2D Flame2D_pState::Strain_Rate(const Flame2D_pState &dWdx,
 					   const Flame2D_pState &dWdy,
 					   const int Axisymmetric,
 					   const Vector2D X){
   Tensor2D strain_rate;
   double radius, div_v;

   /***************** Strain rate (+dilatation) **********************/	
   div_v = dWdx.v.x + dWdy.v.y;
   if (Axisymmetric == AXISYMMETRIC_X) {
     radius = (X.x < MICRO) ? MICRO : X.x;    //fabs(X.x) ??
     div_v += v.x/radius;
   } else if (Axisymmetric == AXISYMMETRIC_Y) {    
     radius = (X.y < MICRO) ? MICRO : X.y;
     div_v += v.y/radius;
   } 

   strain_rate.xx = dWdx.v.x-div_v/THREE;
   strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
   strain_rate.yy = dWdy.v.y-div_v/THREE;

   if (Axisymmetric == PLANAR) {
     strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
   } else if (Axisymmetric == AXISYMMETRIC_X) {
     strain_rate.zz = v.x/radius-div_v/THREE;
   } else if (Axisymmetric == AXISYMMETRIC_Y) {
     strain_rate.zz = v.y/radius-div_v/THREE;
   } 

   return strain_rate;  
 }

 /***************** Laminar (molecular) fluid stress ***********************/
 inline void Flame2D_pState::Laminar_Stress(const Flame2D_pState &dWdx,
 					  const Flame2D_pState &dWdy,
 					  const int Axisymmetric,
 					  const Vector2D X){

   tau = TWO*mu()*Strain_Rate(dWdx, dWdy, Axisymmetric, X); 
 }



//  //Check for unphysical properties
//  /**********************************************************/
//  /* If unphysical properties and using global timestepping */ 
//  /* stop simulation                                        */
//  /**********************************************************/ 

//  inline bool Flame2D_pState::Unphysical_Properties_Check(Flame2D_cState &U, const int Flow_Type, const int n) {
//    if ((Flow_Type == FLOWTYPE_INVISCID ||
//         Flow_Type == FLOWTYPE_LAMINAR) &&
//        (U.rho <= ZERO ||!U.negative_speccheck(n) ||U.es() <= ZERO)) {
//     cout << "\n " << CFFC_Name() 
// 	 << " Flame2D ERROR: Negative Density || Energy || mass fractions: \n"
// 	 << U << "\n";
//     return false;
  
//   }
//   if ((Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
//        Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA) &&
//       (U.rho <= ZERO || !U.negative_speccheck(n) ||U.es() <= ZERO ||
//        U.rhok < ZERO ||U.rhoomega < ZERO)) {
//     cout << "\n " << CFFC_Name() 
// 	 << " Flame2D ERROR: Negative Density || Energy || mass fractions || Turbulent kinetic energy || : \n"
// 	 << " Dissipation rate per unit turbulent kinetic energy || : \n"
// 	 <<U << "\n";
//     return false;
//   }else{
//     return true;
//   }  
// } 
// /**************************************************************
//   Check for -ve mass fractions and set small -ve values
//   to ZERO. Then check that no mass is lost and that 
//   the mass fractions still sum = 1

//   Return "true" if passes
//   and "false" if failed

//   Should add "strict" and "anything goes" flags
//   currently only set to give warnings, but still 
//   continues.

// ***************************************************************/
// inline bool Flame2D_pState::negative_speccheck(void) {
//   double sum(ZERO);
//   double SPEC_TOLERANCE(SPEC_TOLERANCE);

//   //-------- Negative Check ------------//
//   for(int i=0; i<ns-1; i++){
//     if(spec[i].c > ONE){ //check for > 1.0
//       spec[i].c = ONE;
//     } else if(spec[i].c < ZERO){  //check for -ve
//       if(spec[i].c > -SPEC_TOLERANCE){  //check for small -ve and set to ZERO 
// 	spec[i].c = ZERO;
//       } else {
// 	spec[i].c = ZERO;
// 	//#ifdef _DEBUG
// 	cout <<"\n pState -ve mass fraction in "<<specdata[i].Speciesname()<<" "<<
// 	  spec[i].c<<" greater than allowed tolerance of "<<-SPEC_TOLERANCE; 
// 	//#endif
//       }      
//     } 
//     sum += spec[i].c;
//   } 
  
//   //   spec[ns-1].c = (ONE - sum);  //PUSH error into NS-1
//   //Spread error across species
//   spec[ns-1].c = max(ONE- sum, ZERO);
//   sum += spec[ns-1].c;
//   for(int i=0; i<ns; i++){
//     spec[i].c = spec[i].c*(ONE/sum);
//   }

//   return true;
//}

/***** Species Concentrations ******************************/
inline double Flame2D_pState::SpecCon(int i) const{
  //returned in kg/m^3 / kg/mol => mol/m^3
  return (rho)*spec[i].c/(specdata[i].Mol_mass());
}

/***** Species Mole Fractions ******************************/
inline double Flame2D_pState::MoleFrac(int i) const{
  //returned in kmol i / kmol mix
  return spec[i].c*Mass()/(specdata[i].Mol_mass());
}


/******* GIBBS Free Energy ********************************
  Eqn. (10.84) (10.87) Anderson
  Gs = Hs - TS
  Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
  ps = cs(M/Ms)*p

  NEVER USE THIS FUNCTION : If you are looping through the
  species and computing the Gibbs free energy, you do an
  EXPENSIVE T() calculation for every species !!!! Even 
  though T() does not change.
***********************************************************/
inline double Flame2D_pState::Gibbs(int species) const{
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
inline double Flame2D_pState::diedip() const{
  double Rlocal = Rtot();
  return (hprime() - Rlocal)/(rho*Rlocal);
}

inline double Flame2D_pState::diedirho() const{
  double Rlocal = Rtot();
  return -p*(hprime() - Rlocal)/(rho*rho*Rlocal);
}

/**************** Copy *************************************/
inline void Flame2D_pState::Copy(const Flame2D_pState &W){
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  for( int i=0; i<ns; i++) spec[i] = W.spec[i];
  tau = W.tau;
  qflux = W.qflux;
}

//**************** Index Operators *************************/
inline double& Flame2D_pState::operator[](int index) {  
  //  assert( index >= 1 && index <= NUM_VAR_FLAME2D );
  switch(index){  
  case 1:
    return rho;    
  case 2:
    return v.x;
  case 3:
    return v.y;
  case 4:
    return p;
  default :
    return spec[index-NUM_FLAME2D_VAR_SANS_SPECIES-1].c;
  };
}

inline const double& Flame2D_pState::operator[](int index) const {  
  //   assert( index >= 1 && index <= NUM_VAR_FLAME2D );
  switch(index){  
  case 1:
    return rho;    
  case 2:
    return v.x;
  case 3:
    return v.y;
  case 4:
    return p;
  default :
    return spec[index-NUM_FLAME2D_VAR_SANS_SPECIES-1].c;
  };

}
/**************************************************************
  Get max of the min temperature of the lowest region
  and min of the max temperature of the highest region
***************************************************************/
inline void Flame2D_pState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; i++){
    temp = max(specdata[i].Low_range(),temp);
  }
  low_temp_range = temp;  
}

inline void Flame2D_pState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; i++){
    temp = min(specdata[i].High_range(),temp);
  }
  high_temp_range = temp;  
}

/********************************************************
 * Flame2D_pState::U -- Conserved solution state.        *
 ********************************************************/
inline Flame2D_cState Flame2D_pState::U(void) const {
  return U(*this);
}

inline Flame2D_cState Flame2D_pState::U(const Flame2D_pState &W) const{
    Flame2D_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    for(int i=0; i<W.ns; i++){
      Temp.rhospec[i].c = W.rho*W.spec[i].c;
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
    } 
    Temp.E = W.E();
    Temp.tau = W.tau;
    Temp.qflux = W.qflux; 
    return Temp;
}

inline Flame2D_cState U(const Flame2D_pState &W) {
  Flame2D_cState Temp;
  Temp.rho = W.rho;
  Temp.rhov = W.rhov();
  for(int i=0; i<W.ns; i++){
    Temp.rhospec[i].c = W.rho*W.spec[i].c;
    Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
    Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
  }  
  Temp.E = W.E(); 
  Temp.tau = W.tau;
  Temp.qflux = W.qflux; 
  return Temp;
}

/**************************************************************************
********************* FLAME2D_CSTATE CONSTRUCTORS **************************
***************************************************************************/

/***************************************************************
  Set up mass fractions memory and initial values, ie. Species Class  
****************************************************************/
inline void Flame2D_cState::set_initial_values(){
  rhospec_memory();
  for(int i=0; i<ns; i++) rhospec[i].c = rho/ns; 
} 

inline void  Flame2D_cState::set_initial_values(const double &value){
  rhospec_memory();
  for(int i=0; i<ns; i++) rhospec[i].c = value; 
}

//user specified
inline void Flame2D_cState::set_initial_values(double *rhomfrac){
  rhospec_memory();
  for(int i=0; i<ns; i++) rhospec[i].c = rhomfrac[i];
}

//another set using species class
inline void Flame2D_cState::set_initial_values( const Species *rhomfrac){
  rhospec_memory();
  for(int i=0; i<ns; i++) rhospec[i] = rhomfrac[i];
}

/**************************************************************************
  Resize the species array, setting the new values to zero
**************************************************************************/
inline void Flame2D_cState::resize_species() {

  // the new number of species
  int n_new( Flame2D_pState::React.num_species );

  // deallocate first
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES < n_new ) { 
    cerr<<"\n WARNING USING STATIC FLAME2D BUILT WITH "<<STATIC_NUMBER_OF_SPECIES 
	<<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<<n_new<<endl; 
    exit(1); 
  }
#else 
  Deallocate();
#endif  
  
  // set the new value
  ns = n_new;
  
  // reallocate species array and set to zero
  rhospec_memory();
  set_initial_values(ZERO);
     
}

/**************************************************************************
  Resize the species array, keeping the values
**************************************************************************/
inline void Flame2D_cState::resize_species( const int n_old,           // the old numbre of species
					   const string* species_old) // the old species names array
{

  // temporary storage
  double* c_old = new double[n_old];
  int index;

  // the new number of species
  int n_new( Flame2D_pState::React.num_species );

  // store the old values
  for (int k=0; k<n_old; k++) c_old[k] = rhospec[k].c;

  // deallocate
#ifdef STATIC_NUMBER_OF_SPECIES
  if( STATIC_NUMBER_OF_SPECIES < n_new ) { 
    cerr<<"\n WARNING USING STATIC FLAME2D BUILT WITH "<<STATIC_NUMBER_OF_SPECIES 
	<<" SPECIES PREDEFINED, HOWEVER ASKING FOR "<<n_new<<endl; 
    exit(1); 
  }
#else 
  Deallocate();
#endif  
  
  // set the new value
  ns = n_new;
  
  // reallocate species array and set to zero
  rhospec_memory();
  set_initial_values(ZERO);
    
  // match up the species and set their respective values
  for (int k=0; k<n_old; k++) {
    index = Flame2D_pState::React.SpeciesIndex(species_old[k]);
    if (index>=0) rhospec[index].c = c_old[k];
    else {
      cerr << "\n Flame2D_cState::resize_species() - Species " 
	   << species_old[k]
	   << " not found in current reaction mechanism.\n";
      exit(-1);
    }
  }// endfor - species
  
  // delete temporary memory
  delete[] c_old;

}


/**************** Copy *************************************/
inline void Flame2D_cState::Copy(const Flame2D_cState &U){
  rho = U.rho;
  rhov = U.rhov; 
  E = U.E; 

  for( int i=0; i<ns; i++){ 
    rhospec[i] = U.rhospec[i];
  } 
  tau = U.tau;
  qflux = U.qflux; 
}

 /**********************************************************************
  * Flame2D_cState::set_gravity -- Set the acceleration due to gravity  *
  *                               in m/s^2.  It acts downwards in the  *
  *                               z-dir (g <= 0)                       *
  **********************************************************************/
inline void Flame2D_cState::set_gravity(const double &g) { // [m/s^2]

  // if gravity is acting upwards (wrong way)
  if (g>0) {
    cerr<<"\n Flame2D_pState::set_gravity() - Gravity acting upwards!!!! \n";
    exit(1);
    
  // gravity acting downwards (-ve), OK
  } else {
    gravity_z = g;
  }
}


/**************** Velocity *********************************/
inline Vector2D Flame2D_cState::v() const{
  return (rhov/rho);
}

/**************** Pressure *********************************/
inline double Flame2D_cState::p() const{
  return (rho*Rtot()*T());
}

//----------------- Index Operator ---------------------//
inline double& Flame2D_cState::operator[](int index) {
  //  assert( index >= 1 && index <= NUM_VAR_FLAME2D );
  switch(index){  
  case 1:
    return rho;    
  case 2:
    return (rhov.x);
  case 3:
    return (rhov.y);
  case 4:
    return (E);
  default :
    return rhospec[index-NUM_FLAME2D_VAR_SANS_SPECIES-1].c;
  };
}

inline const double& Flame2D_cState::operator[](int index) const{
  //  assert( index >= 1 && index <= NUM_VAR_FLAME2D ); 
  switch(index){  
  case 1:
    return rho;    
  case 2:
    return (rhov.x);
  case 3:
    return (rhov.y);
  case 4:
    return (E);
  default :
    return rhospec[index-NUM_FLAME2D_VAR_SANS_SPECIES-1].c;
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
inline bool Flame2D_cState::negative_speccheck(const int &step) {
  double sum(ZERO);
  double temp(ZERO);

  //-------- Negative Check ------------//     
  for(int i=0; i<ns-1; i++){
    temp = rhospec[i].c/rho;    
    if(temp > ONE){            //check for > 1.0
      rhospec[i].c =rho;
      temp = ONE; 

    } else if(temp < ZERO){   //check for -ve
      if(temp > -SPEC_TOLERANCE){  //check for small -ve and set to ZERO 
	rhospec[i].c = ZERO;
	temp = ZERO;
      } else {

#ifdef _DEBUG
	cout<<"\ncState -ve mass fraction in "<<specdata[i].Speciesname()<<" "<<
	  temp<<" greater than allowed tolerance of "<<-SPEC_TOLERANCE; 
#endif

	if( step < 10){
	  return false;
	} else { 

#ifdef _DEBUG
	  cout<<"\ncState rhospec["<<i<<"] = "<<rhospec[i].c<<" -ve mass fraction larger than tolerance,"
	      <<" but setting to zero and continuing anyway. ";
#endif
	  rhospec[i].c = ZERO;
	  temp = ZERO;	 
	}
      }
    } else {
      // cout<<"\n negative_speccheck else";
    }
    sum += temp;
  } 

//   rhospec[ns-1].c = rho*(ONE - sum); //PUSH error into NS-1 species, probably N2

  temp = max(ONE- sum, ZERO);    //Spread Error across species
  sum += temp;
  rhospec[ns-1].c = rho*temp;
  for(int i=0; i<ns; i++){
    rhospec[i].c = rhospec[i].c*(ONE/sum);
  } 
  
  return true;
}

// USED IN MULTIGRID
inline bool Flame2D_cState::Unphysical_Properties(){
  Flame2D_cState::Unphysical_Properties_Check(10);
}

/**************************************************************
  Unphysical_Properties_Check
 ***************************************************************/
inline bool Flame2D_cState::Unphysical_Properties_Check(const int n){

  // check for nan's, inf's etc.... debugging !!!
#ifdef _DEBUG
  for( int i=0; i<NUM_FLAME2D_VAR_SANS_SPECIES+ns; i++){ 
    if( this[i] != this[i]){ cout<<"\n nan's in solution, variable "<<i<<endl; exit(1); return false; }
  }
#endif

  if (rho <= ZERO ||!negative_speccheck(n) ||es() <= ZERO) {
    cout << "\n " << CFFC_Name() 
	 << " Flame2D ERROR: Negative Density || Energy || mass fractions: \n" << *this <<endl;
    return false;
  }else{
    return true ;
  }
} 

/**************************************************************
  Sum N-1 species.
***************************************************************/
inline double Flame2D_cState::sum_species(void) const{
  double sum=ZERO;
  for(int i=0; i<ns-1; i++){
    sum += rhospec[i].c;
  }
  return sum/rho;
}

/**********************************************************************
 * Flame2D_cState::Copy_Multigrid_State_Variables --                   *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void Flame2D_cState::Copy_Multigrid_State_Variables(const Flame2D_cState &Ufine) {
  Copy(Ufine);
  Zero_Non_Multigrid_State_Variables();
} 
  
/**********************************************************************
 * Flame2D_cState::Zero_Non_Multigrid_State_Variables --               *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void Flame2D_cState::Zero_Non_Multigrid_State_Variables(void) {
}

/**************************************************************
  Get max and min temperature ranges for data
***************************************************************/
inline void Flame2D_cState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; i++){
    temp = max(specdata[i].Low_range(),temp);
  } 
  low_temp_range = temp;  
}

inline void Flame2D_cState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; i++){
    temp = min(specdata[i].High_range(),temp);
  } 
  high_temp_range = temp;  
}


/********************************************************
 * Flame2D_cState::W -- Primitive solution state.       *
 ********************************************************/
inline Flame2D_pState Flame2D_cState::W(void) const {
  return W(*this);
}

inline Flame2D_pState Flame2D_cState::W(const Flame2D_cState &U) const{
    Flame2D_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    for(int i=0; i<U.ns; i++){
      Temp.spec[i].c = U.rhospec[i].c/U.rho;
      Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
      Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
    }   
    Temp.p = U.p();
    Temp.tau = U.tau;
    Temp.qflux = U.qflux; 

    return Temp;
}

inline Flame2D_pState W(const Flame2D_cState &U) {
  Flame2D_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  for(int i=0; i<U.ns; i++){
    Temp.spec[i].c = U.rhospec[i].c/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  } 
  Temp.p = U.p();
  Temp.tau = U.tau;
  Temp.qflux = U.qflux;

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
extern Flame2D_pState Reflect(const Flame2D_pState &W,
	      	              const Vector2D &norm_dir);

extern Flame2D_pState Free_Slip(const Flame2D_pState &Win,
			     const Flame2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG);

extern Flame2D_pState No_Slip(const Flame2D_pState &Win,
			     const Flame2D_pState &Wout,
			     const Vector2D &norm_dir,
			     const int &TEMPERATURE_BC_FLAG);

extern Flame2D_pState Moving_Wall(const Flame2D_pState &Win,
				 const Flame2D_pState &Wout,
				 const Vector2D &norm_dir,				 
				 const double &wall_velocity,
				 const int &TEMPERATURE_BC_FLAG);

extern Flame2D_pState BC_Characteristic_Pressure(const Flame2D_pState &Wi,
						const Flame2D_pState &Wo,
						const Vector2D &norm_dir);

extern Flame2D_pState BC_1DFlame_Inflow(const Flame2D_pState &Wi,
				     const Flame2D_pState &Wo, 
				     const Flame2D_pState &Woutlet,
				     const Vector2D &norm_dir);

extern Flame2D_pState BC_2DFlame_Inflow(const Flame2D_pState &Wi,
				       const Flame2D_pState &Wo, 				
				       const Vector2D &norm_dir);

extern Flame2D_pState BC_1DFlame_Outflow(const Flame2D_pState &Wi,
					const Flame2D_pState &Wo,
					const Flame2D_pState &Winlet,
					const Vector2D &norm_dir);

extern Flame2D_pState BC_2DFlame_Outflow(const Flame2D_pState &Wi, 
					const Flame2D_pState &Wo,
					const Vector2D &norm_dir);

/*******************************************************
 * Exact Test Case Solution Functions                  *
 *******************************************************/
extern Flame2D_pState RinglebFlow(const Flame2D_pState &Wdum,
				  const Vector2D X);

extern Flame2D_pState ViscousChannelFlow(const Flame2D_pState &Wdum,
					const Vector2D X,
					const double Vwall,
					const double dp);

extern Flame2D_pState FlatPlate(const Flame2D_pState &Winf,
			       const Vector2D X,
			       double &eta,
			       double &f,
			       double &fp,
			       double &fpp);

/*******************************************************
 * External Flux Function Functions                    *
 *******************************************************/
extern Flame2D_pState RoeAverage(const Flame2D_pState &Wl,
	      	                 const Flame2D_pState &Wr);

// HLLE
extern Flame2D_cState FluxHLLE_x(const Flame2D_pState &Wl,
				const Flame2D_pState &Wr,
				const int &Preconditioning);

extern Flame2D_cState FluxHLLE_x(const Flame2D_cState &Ul,
				const Flame2D_cState &Ur,
				const int &Preconditioning);
  
extern Flame2D_cState FluxHLLE_n(const Flame2D_pState &Wl,
				const Flame2D_pState &Wr,
				const Vector2D &norm_dir,
				const int &Preconditioning);

extern Flame2D_cState FluxHLLE_n(const Flame2D_cState &Ul,
				const Flame2D_cState &Ur,
				const Vector2D &norm_dir,
				const int &Preconditioning);

// Linde
extern Flame2D_cState FluxLinde(const Flame2D_pState &Wl,
			       const Flame2D_pState &Wr);

extern Flame2D_cState FluxLinde(const Flame2D_cState &Ul,
			       const Flame2D_cState &Ur);


extern Flame2D_cState FluxLinde_n(const Flame2D_pState &Wl,
				 const Flame2D_pState &Wr,
				 const Vector2D &norm_dir);

extern Flame2D_cState FluxLinde_n(const Flame2D_cState &Ul,
				 const Flame2D_cState &Ur,
				 const Vector2D &norm_dir);

// Roe
extern Flame2D_pState WaveSpeedPos(const Flame2D_pState &lambda_a,
				  const Flame2D_pState &lambda_l,
				  const Flame2D_pState &lambda_r);
			

extern Flame2D_pState WaveSpeedNeg(const Flame2D_pState &lambda_a,
				  const Flame2D_pState &lambda_l,
				  const Flame2D_pState &lambda_r);

extern Flame2D_pState WaveSpeedAbs(const Flame2D_pState &lambda_a,
				  const Flame2D_pState &lambda_l,
				  const Flame2D_pState &lambda_r);
 
extern Flame2D_pState HartenFixPos(const Flame2D_pState &lambda_a,
				  const Flame2D_pState &lambda_l,
				  const Flame2D_pState &lambda_r );
				

extern Flame2D_pState HartenFixNeg(const Flame2D_pState &lambda_a,
				  const Flame2D_pState &lambda_l,
				  const Flame2D_pState &lambda_r);
				  

extern Flame2D_pState HartenFixAbs(const Flame2D_pState &lambdas_a,
				  const Flame2D_pState &lambdas_l,
				  const Flame2D_pState &lambdas_r );
				

extern Flame2D_cState FluxRoe_x(const Flame2D_pState &Wl,
			       const Flame2D_pState &Wr,
			       const int &Preconditioning, 
			       const int &flow_type_flag,
			       const double &deltax);

extern Flame2D_cState FluxRoe_x(const Flame2D_cState &Ul,
			       const Flame2D_cState &Ur,
			       const int &Preconditioning, 
			       const int &flow_type_flag,
			       const double &deltax);

extern Flame2D_cState FluxRoe_n(const Flame2D_pState &Wl,
			       const Flame2D_pState &Wr,
			       const Vector2D &norm_dir,
			       const int &Preconditioning, 
			       const int &flow_type_flag,
			       const double &deltax);

extern Flame2D_cState FluxRoe_n(const Flame2D_cState &Ul,
			       const Flame2D_cState &Ur,
			       const Vector2D &norm_dir,
			       const int &Preconditioning, 
			       const int &flow_type_flag,
			       const double &deltax);

extern Flame2D_cState FluxAUSMplus_up(const Flame2D_pState &Wl,
				     const Flame2D_pState &Wr);

extern Flame2D_cState FluxAUSMplus_up(const Flame2D_cState &Wl,
				     const Flame2D_cState &Wr);

extern Flame2D_cState FluxAUSMplus_up_n(const Flame2D_pState &Wl,
				       const Flame2D_pState &Wr,
				       const Vector2D &norm_dir);

extern Flame2D_cState FluxAUSMplus_up_n(const Flame2D_cState &Wl,
				       const Flame2D_cState &Wr,
				       const Vector2D &norm_dir);

/* Viscous Solution flux (laminar+turbulent) */
extern Flame2D_cState Viscous_FluxArithmetic_n(const Flame2D_cState &Ul,
   			                      const Flame2D_pState &dWdx_l,
			                      const Flame2D_pState &dWdy_l,
				              const Flame2D_cState &Ur,
   			                      const Flame2D_pState &dWdx_r,
			                      const Flame2D_pState &dWdy_r,
				              const Vector2D &norm_dir);

extern Flame2D_cState Viscous_Flux_n(Flame2D_pState &W,
				    const Flame2D_pState &dWdx,
				    const Flame2D_pState &dWdy,
				    const int Axisymmetric,
				    const Vector2D X,			
				    const Vector2D &norm_dir);

extern double WallShearStress(const Flame2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &norm_dir);

/*********************************************************/
extern Flame2D_pState Rotate(const Flame2D_pState &W,
	      	             const Vector2D &norm_dir);


extern Vector2D HLLE_wavespeeds(const Flame2D_pState &Wl,
                                const Flame2D_pState &Wr,
                                const Vector2D &norm_dir); 

#endif //end _FLAME2D_STATE_INCLUDED 
