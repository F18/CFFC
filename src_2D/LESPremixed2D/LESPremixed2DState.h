/****************** LESPremixed2DState.h **************************************
  This class defines the state variables and constructors for the 
  2D Favre-Filtered Navier-Stokes with Multiple Species solution.

  NOTES:  
         - be careful with set_species_data(), set_initial_values() and
           set_initial_values_scal() due to dynamic memory 

         - if #define STATIC_LESPREMIXED2D_SPECIES is set dynamic memory
           is not used so code is faster, however code is not as flexible
           as only up to STATIC_LESPREMIXED2D_SPECIES can be used without recompling.
  
         - Avoid cState.T() as much as possible as it is $$$ to calculate, pass
           the temperature from pState if possible as it is cheap (p/rhoR)             
              
         - The Static "specdata" has to be deleted outside the class
           in order to avoid seg faults. In LESPremixed2D it is done by calling
           the deallocate function from the Input destructor as it is kept
           until the end of the program.

         - if #define THICKENED_FLAME_ON is set the thickened flame and power-law
           stuff will be activated

  TODO:  - Unphysical_properties & negativespec_check functions need to be  
           worked through

                   
***********************************************************************/
#ifndef _LESPREMIXED2D_STATE_INCLUDED 
#define _LESPREMIXED2D_STATE_INCLUDED

// define early so can be used in other classes
class LESPremixed2D_cState;
class LESPremixed2D_pState;

// Required C++ libraries
#include <iostream>

using namespace std;

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
#include "../Physics/Species.h"
#endif //_SPECIES_INCLUDED

#ifndef _NASARP1311_DATA_INCLUDED
#include "../Physics/NASAData/NASARP1311data.h"
#endif


// Other header files
#ifndef _REACTIONS_INCLUDED
#include "../Reactions/Reactions.h"
#endif //_REACTIONS_INCLUDED

#ifndef _POWER_LAW_INCLUDED
#include "PowerLaw.h"
#endif //_POWER_LAW_INCLUDED

#ifndef _SCALARS_INCLUDED
#include "Scalars.h"
#endif //_SCALARS_INCLUDED

#ifndef _SFS_MODELLING_INCLUDED
#include "../Turbulent2D/SFSModelling.h"
#endif // _SFS_MODELLING_INCLUDED



//Temperature convergence tolerance in
//LESPremixed2D_cState::T(void)
// these should be moved to CFD.h or Math.h
#define CONV_TOLERANCE 1e-8  //Used for temperature convergence
#define gravity_z -9.81      //m/s^2 acceleration due to gravity  
#define SPEC_TOLERANCE 1e-8  //Used in negative_speccheck for species round off (was MICRO)      

//number of fixed variables in the LESPremixed2D class
#define NUM_LESPREMIXED2D_VAR_SANS_SPECIES 4  //rho, v(2), p


// If you define this variable, the number of species will be
// predetermined for faster calculations.., however it is not as general 
//#define STATIC_LESPREMIXED2D_SPECIES 6 //2 AIR, 6 2STEP_CH4
#define STATIC_LESPREMIXED2D_SPECIES 5 //2 AIR, 6 2STEP_CH4 5 1STEP_CH4


// If you define this variable, the thickened flame model will be on 
//#define THICKENED_FLAME_ON



/*!
 * Class: LESPremixed2D_pState
 *
 * @brief  Primitive variable solution state class definition for an
 *         inviscid, laminar, or turbulent premixed reacting gas-flow.
 *
 * Primitive variable solution state class definition for an inviscid, 
 * laminar, or turbulent premixed reacting gas-flow.
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
class LESPremixed2D_pState {

   private:
   //all public ....
   protected:
   public:
    //@{ @name Primitive variables and associated constants:
   double         rho;   //!< Density.
   Vector2D         v;   //!< Flow velocity (2D)  
   double           p;   //!< Pressure.
#ifdef STATIC_LESPREMIXED2D_SPECIES
   Species       spec[STATIC_LESPREMIXED2D_SPECIES];
#else 
   Species      *spec;   //!< Species class c[ns]
#endif
   double     *scalar;   //!< Scalars using scalar[nscal]     

   Tensor2D                    tau; //!< Shear Stress Tensor
   Vector2D                  qflux; //!< Heat Flux Vector  
   Tensor2D                 lambda; //!< SFS Stress Tensor
   Vector2D                  theta; //!< Turbulent Heat Flux Vector  

#ifdef THICKENED_FLAME_ON
   PowerLaw                  flame; //!< SFS wrinkling factor and thickening factor
#endif 

  //! Static Variables 
  static int                         ns; //!< number of species
  static int                      nscal; //!< number of scalars
  static NASARP1311data       *specdata; //!< Global Species Data
  static double                *Schmidt; //!< Schmidt Number for each species
  static Reaction_set             React; //!< Global Reaction Data
  static Set_scalar            Scal_sys; //!< Set the group of scalars to be solved in the model
  static double          low_temp_range; //!< Low temp data range
  static double         high_temp_range; //!< High temp data range
  static int      NUM_VAR_LESPREMIXED2D; //!< Number of LESPremixed2D variables (4+ns)
  static double                    Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
  static SubfilterScaleModels  SFSmodel; //!< Subfilter scale model
  static double            filter_width; //!< Constant filter width
  //@}

  
  //@{ @name Premixed flame parameters:
  static double      laminar_speed; //!< Unstrained laminar flame speed 
  static double  laminar_thickness; //!< Unstrained laminar flame thickness
  static double            TFactor; //!< Maximum thickening factor
  static double     adiabatic_temp; //!< Adiabatic temperature
  static double  equivalence_ratio; //!< Equivalence ratio
  static double      reactants_den; //!< Reactants density
  //@}
  

  
  //default constructors of many flavours, hopefully one is right 4U   
  //v.zero(); tau.zero(); qflux.zero(); lambda.zero(); theta.zero();  //Unecessary as Vector2D and Tensor2D defaults are ZERO
  //@{ @name Creation, copy, and assignment constructors.
 
   LESPremixed2D_pState(): rho(DENSITY_STDATM), p(PRESSURE_STDATM)  
                 { scalar = NULL; set_initial_values_scal();  specnull();  set_initial_values(); }

   LESPremixed2D_pState(const double &value): rho(value), v(value), p(value)
                { scalar = NULL; set_initial_values_scal();  specnull(); set_initial_values(value); }

   LESPremixed2D_pState(const double &d, const Vector2D &V, const double &pre):
                rho(d), v(V), p(pre) 
 		{ scalar = NULL; set_initial_values_scal();  specnull();  set_initial_values(); }

   LESPremixed2D_pState(const double &d, const double &vx, const double &vy, const double &pre):
		rho(d), v(vx,vy), p(pre) 
 		{ scalar = NULL; set_initial_values_scal();  specnull();  set_initial_values(); }

   LESPremixed2D_pState(const double &d, const Vector2D &V, const double &pre, const double *scal):
		rho(d), v(V), p(pre) 
     { scalar = NULL; set_initial_values_scal(scal);  specnull();  set_initial_values(); } 

   LESPremixed2D_pState(const double &d, const Vector2D &V, const double &pre, const double &value):
		rho(d), v(V), p(pre) 
 		{ scalar = NULL; set_initial_values_scal();  specnull();  set_initial_values(value); } 

   LESPremixed2D_pState(const double &d, const double &vx, const double &vy, const double &pre, const double *scal):
                 rho(d), v(vx,vy), p(pre) 
     { scalar = NULL; set_initial_values_scal(scal);  specnull();  set_initial_values(); }

   LESPremixed2D_pState(const double &d, const double &vx, const double &vy, const double &pre, const double &value):
                 rho(d), v(vx,vy), p(pre)                
 		{ scalar = NULL; set_initial_values_scal();  specnull();  set_initial_values(value); }

   LESPremixed2D_pState(const double &d, const double &vx, const double &vy, const double &pre, 
                 const double *scal, const double &value):
                 rho(d), v(vx,vy), p(pre) 
 		{ scalar = NULL; set_initial_values_scal(scal);  specnull();  set_initial_values(value); }

   LESPremixed2D_pState(const double &d, const double &vx, const double &vy, const double &pre, 
                 const double *scal, const Species *mfrac): 
                 rho(d), v(vx,vy), p(pre) 
 		{ scalar = NULL; set_initial_values_scal(scal);  specnull();  set_initial_values(mfrac); }

   LESPremixed2D_pState(const double &d, const Vector2D &V, const double &pre,
                 const double *scal, const Species *mfrac): 
                 rho(d), v(V), p(pre) 
                 { scalar = NULL; set_initial_values_scal(scal);  specnull(); set_initial_values(mfrac); }

   //this is needed for the operator overload returns!!!!
   LESPremixed2D_pState(const LESPremixed2D_pState &W): rho(W.rho), v(W.v), p(W.p),
 				 tau(W.tau), qflux(W.qflux), lambda(W.lambda), theta(W.theta) 
#ifdef THICKENED_FLAME_ON
		   ,flame(W.flame)
#endif 
		   { scalar = NULL; set_initial_values_scal(W.scalar); 
		     specnull(); set_initial_values(W.spec); }            
  //@}

   //read in ns species data, call only once as its static
   void set_species_data(const int &, const int &, const string *, const char *,
 			 const double&, const double *, const int &);

   //set initial data values predominately used internally !!!
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *cfrac);
   void set_initial_values(const Species *mfrac);

   void set_initial_values_scal();
   void set_initial_values_scal(const double &value);
   void set_initial_values_scal(const double *scal);


   //Copy constructor, cheaper than = operator
   void Copy(const LESPremixed2D_pState &W);

   // return the number of variables - number of species
   int NumVarSansSpecies() const { return NUM_VAR_LESPREMIXED2D-ns; }

   /*************** VACUUM OPERATOR *********************/
   void Vacuum(){ rho=ZERO; v.zero(); p=ZERO;
     if (nscal) for(int i=0; i<nscal; ++i) scalar[i] = ZERO;
     for(int i=0; i<ns; ++i)  spec[i].Vacuum();
     tau.zero(); qflux.zero(); lambda.zero(); theta.zero(); 
   }

   void zero_non_sol(){
     for(int i=0; i<ns; ++i){
       spec[i].gradc.zero();
       spec[i].diffusion_coef=ZERO;
     }
     tau.zero(); qflux.zero(); lambda.zero(); theta.zero(); 
   }  
 

  //! Set premixed flame static variables.
  void set_premixed_flame_variables(const double &lam_thickness,
				    const double &lam_speed,
				    const double &thickening_factor,
                                    const double &adia_temp,
                                    const double &equi_ratio,
                                    const double &reac_den);

  //! Set SFS modelling variables.
  void set_SFSmodel_variables(const double &delta,
			      const double &Smagorinsky_coef,
			      const double &Yoshizawa_coef);
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
   double pmodified(void) const;   //Turbulence modified pressure
   double gamma_guess(void) const; //gamma for iterative schemes
   double a(void);                 //speed of sound
   double a(void) const;
   double amodified(void) const;   //Turbulence modified speed of sound
   double k(void) const;           //SFS turbulence kinetic energy
   bool negative_speccheck(void) ; //-ve mass frac check and sets small -ve c's to ZERO  
//    bool Unphysical_Properties_Check(LESPremixed2D_cState &U, const int Flow_Type, const int n) ;   
   double SpecCon(int i) const;     //Species i concentration (rho*c/mol_mass)
   double Gibbs(int species) const; //Gibbs Free Energy (H-TS) for species

   /**************** turbulence model related parameters*********/
   double mu_t(const Tensor2D &strain_rate, const int &Flow_Type) const;      
   double Pr_turb(void) const;      
   double Sc_turb(void) const;      
   double Kappa_turb(const double &mut) const;      
   double Dm_turb(const double &mut) const;

   /************* Dimensionless Parameters ********************/
   double Schmidt_No(const int &) const;
   double Prandtl() const;
   double Lewis(const int &) const;

   /************** Temperature Derivatives *******************/
   double diedip() const;
   double diedirho() const;
   double dmudT(void) const;
   double dkappadT(void) const;

   /************ Strain rate tensor, laminar stress and SFS stress tensors ***********/
   Tensor2D Strain_Rate(const LESPremixed2D_pState &dWdx,
			const LESPremixed2D_pState &dWdy,
			const int Flow_Type,
			const int Axisymmetric,
			const Vector2D X);

   void Laminar_Stress(const Tensor2D &strain_rate);

   void SFS_Stress(const Tensor2D &strain_rate, const int &Flow_Type);

   /************ Heat Flux vector thermal Diffusion ***********/
   Vector2D thermal_diffusion(void) const;

   /*************** Conserved solution state. ****************/
   LESPremixed2D_cState U(void) const;
   LESPremixed2D_cState U(const LESPremixed2D_pState &W) const;
   friend LESPremixed2D_cState U(const LESPremixed2D_pState &W);

   /**************** Fluxes ***********************************/
   /* Inviscid Solution flux (x-direction). */
   LESPremixed2D_cState Fx(void) const;

   /*********** Inviscid Flux Jacobian X ***********************/
   friend void dFIdU(DenseMatrix &dFdU, const LESPremixed2D_pState &W, const int Flow_Type);
   friend void dFIdU_FD(DenseMatrix &dFdU, const LESPremixed2D_pState &W, const int Flow_Type);
   friend void dFIdW(DenseMatrix &dFdW, const LESPremixed2D_pState &W, const int Flow_Type);
   friend void dFIdW_FD(DenseMatrix &dFdW, const LESPremixed2D_pState &W, const int Flow_Type);

   /** Conserved/Primitive  & Primitive/Conservered Jacobians **/
   void dWdU(DenseMatrix &dWdQ, const int Flow_Type)const;
   void dWdU_FD(DenseMatrix &dWdQ, const int Flow_Type);

   void dUdW(DenseMatrix &dQdW, const int Flow_Type);
   void dUdW_FD(DenseMatrix &dUdW, const int Flow_Type);

   /************* Eigenvalues *********************************/
   /* Eigenvalue(s) (x-direction). */ 
   LESPremixed2D_pState lambda_x(void) const;

   /************** Precondtioned Eigenvalues ******************/
   LESPremixed2D_pState lambda_preconditioned_x(const double &MR2) const;

   /**************** Eigenvectors *****************************/  
   LESPremixed2D_cState rc_x(const int &index) const; // Conserved right (x-direction)
   LESPremixed2D_pState lp_x(const int &index) const; // Primitive left (x-direction)

   /************** Preconditioned Eigenvectors ****************/
   LESPremixed2D_cState rc_x_precon(const int &index,const double &MR2) const;  // Conserved right (x-direction)
   LESPremixed2D_pState lp_x_precon(const int &index,const double &MR2) const;  // Primitive left  (x-direction)

   /*************** Preconditioner ****************************/
   double u_plus_aprecon(const double &u,const int &flow_type_flag,const double &deltax) const;
   void u_a_precon(const double &UR,double &uprimed, double &cprimed) const;
   double Mr2(const int &flow_type_flag, const double &deltax) const;
   void Low_Mach_Number_Preconditioner(DenseMatrix &P,const int &flow_type_flag, const double &deltax) const; 
   void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,const int &flow_type_flag, const double &deltax) const; 

   /********** Axisymmetric Source Terms **********************/
   LESPremixed2D_cState Sa_inviscid(const Vector2D &X,
				    const int Flow_Type,
				    const int Axisymmetric) const;

   LESPremixed2D_cState Sa_viscous(const LESPremixed2D_pState &dWdx, 
				   const LESPremixed2D_pState &dWdy, 
				   const Vector2D &X,
				   const int Flow_Type,
				   const int Axisymmetric);

   //Due to axisymmetric coordinate system
  void dSa_idU(DenseMatrix &dSa_idU, const Vector2D &X, const int Flow_Type,const int Axisymmetric) const;  //Inviscid source Jacobian
  void dSa_idU_FD(DenseMatrix &dSa_idU, const Vector2D &X, const int Flow_Type,const int Axisymmetric) const;
  void dSa_vdW(DenseMatrix &dSa_VdW, const LESPremixed2D_pState &dWdx,
 	       const LESPremixed2D_pState &dWdy,const Vector2D &X, 
 	       const int Flow_Type,const int Axisymmetric, 
 	       const double d_dWdx_dW, const double d_dWdy_dW) const;  //Viscous source Jacobian

   /****** Source terms associated with finite-rate chemistry ******/
   LESPremixed2D_cState Sw(int &REACT_SET_FLAG, const int Flow_Type ) const;
   void dSwdU(DenseMatrix &dSwdU, const int &Flow_Type,const int &solver_type) const; //Jacobian
   void dSwdU_FD(DenseMatrix &dSwdU, const int Flow_Type) const;
   double dSwdU_max_diagonal(const int &Preconditioned,					 
			     const int &Viscous,
			     const double &delta_n,
			     const int &solver_type) const;

   /****** Source terms associated with gravitational forces ******/
   LESPremixed2D_cState Sg(void) const;
   void dSgdU(DenseMatrix &dSgdU) const; //Jacobian

   /********** Source terms associated with turbulence model *****************/
   LESPremixed2D_cState S_turbulence_model(const LESPremixed2D_pState &dWdx, 
					   const LESPremixed2D_pState &dWdy, 
					   const Vector2D &X,
					   const int Flow_Type,
					   const int Axisymmetric);

  /****** Source terms associated with dual time stepping *******/
  LESPremixed2D_cState S_dual_time_stepping(const LESPremixed2D_cState &U,
					    const LESPremixed2D_cState &Ut,
					    const LESPremixed2D_cState &Uold,
					    const double &dTime,
					    const int &first_step) const;

  /************** Premixed combustion ************************/
  LESPremixed2D_pState premixed_mfrac(const LESPremixed2D_pState &Wo);
  int FlameJumpLowMach_x(LESPremixed2D_pState &Wu,
			 LESPremixed2D_pState &Wb,
			 const int &Flow_Type);
  int FlameJump_x(const LESPremixed2D_pState &Wu,
		  LESPremixed2D_pState &Wb,
		  const int &Flow_Type);
  int FlameJump_n(const LESPremixed2D_pState &Wu,
		  LESPremixed2D_pState &Wb,
		  const Vector2D &norm_dir,
		  const int &Flow_Type);

  /************** FSD Model Source Term **********************/
  double HeatRelease_Parameter(void)const;
  double SFS_Kinetic_Energy_Fsd(const LESPremixed2D_pState &dWdx,
                                const LESPremixed2D_pState &dWdy,
                                const int &Flow_Type,
                                const Tensor2D &strain_rate)const;
  double Efficiency_Function_Fsd(const LESPremixed2D_pState &dWdx,
                                 const LESPremixed2D_pState &dWdy,
                                 const int &Flow_Type,
                                 const Tensor2D &strain_rate) const; 
  double Progvar_Species_Grad(void) const;
  double Reaction_Rate_Progvar(const LESPremixed2D_pState &dWdx,
                               const LESPremixed2D_pState &dWdy) const;
  double Reaction_Rate_Algebraic(const LESPremixed2D_pState &dWdx,
                                 const LESPremixed2D_pState &dWdy,
                                 const int &Flow_Type,
                                 const Tensor2D &strain_rate) const;
  double Reaction_Rate_NGT_C_Fsd(const LESPremixed2D_pState &dWdx,
                                 const LESPremixed2D_pState &dWdy) const;
  double Reaction_Rate_Fsd(const LESPremixed2D_pState &dWdx,
                           const LESPremixed2D_pState &dWdy) const;
  double M_x(const LESPremixed2D_pState &dWdx,
             const LESPremixed2D_pState &dWdy) const;
  double M_y(const LESPremixed2D_pState &dWdx,
             const LESPremixed2D_pState &dWdy) const;
  double Resolved_Strain(const LESPremixed2D_pState &dWdx,
                         const LESPremixed2D_pState &dWdy) const;
  double Resolved_Propagation_Curvature(const LESPremixed2D_pState &dWdx,
                                        const LESPremixed2D_pState &dWdy) const;
  double SFS_Strain(const LESPremixed2D_pState &dWdx,
                    const LESPremixed2D_pState &dWdy,
                    const int &Flow_Type,
                    const Tensor2D &strain_rate) const;
  double SFS_Curvature(const LESPremixed2D_pState &dWdx,
                       const LESPremixed2D_pState &dWdy,
                       const int &Flow_Type) const;
  double M_xx(const LESPremixed2D_pState &dWdx,
              const LESPremixed2D_pState &dWdy,
              const LESPremixed2D_pState &d_dWdx_dx,
              const LESPremixed2D_pState &d_dWdx_dy,
	      const LESPremixed2D_pState &d_dWdy_dy) const;
  double M_xy(const LESPremixed2D_pState &dWdx,
              const LESPremixed2D_pState &dWdy,
              const LESPremixed2D_pState &d_dWdx_dx,
              const LESPremixed2D_pState &d_dWdx_dy,
	      const LESPremixed2D_pState &d_dWdy_dy) const;
  double M_yy(const LESPremixed2D_pState &dWdx,
              const LESPremixed2D_pState &dWdy,
              const LESPremixed2D_pState &d_dWdx_dx,
              const LESPremixed2D_pState &d_dWdx_dy,
	      const LESPremixed2D_pState &d_dWdy_dy) const;
  double Resolved_Curvature(const LESPremixed2D_pState &dWdx,
                            const LESPremixed2D_pState &dWdy,
                            const LESPremixed2D_pState &d_dWdx_dx,
                            const LESPremixed2D_pState &d_dWdx_dy,
                            const LESPremixed2D_pState &d_dWdy_dy) const;
  double Resolved_Propagation(const LESPremixed2D_pState &dWdx,
                              const LESPremixed2D_pState &dWdy,
                              const LESPremixed2D_pState &d_dWdx_dx,
                              const LESPremixed2D_pState &d_dWdx_dy,
                              const LESPremixed2D_pState &d_dWdy_dy) const;
  double Resolved_Convection_Progvar (const LESPremixed2D_pState &dWdx,
                                      const LESPremixed2D_pState &dWdy) const;
  double Resolved_Convection_Fsd (const LESPremixed2D_pState &dWdx,
                                  const LESPremixed2D_pState &dWdy) const;
  double NGT_Progvar (const LESPremixed2D_pState &dWdx,
                      const LESPremixed2D_pState &dWdy) const;
  double NGT_Fsd (const LESPremixed2D_pState &dWdx,
                  const LESPremixed2D_pState &dWdy,
                  const LESPremixed2D_pState &d_dWdx_dx,
                  const LESPremixed2D_pState &d_dWdx_dy,
                  const LESPremixed2D_pState &d_dWdy_dy) const;
  double Heat_Release_Strain (const LESPremixed2D_pState &dWdx,
                              const LESPremixed2D_pState &dWdy,
                              const LESPremixed2D_pState &d_dWdx_dx,
                              const LESPremixed2D_pState &d_dWdx_dy,
                              const LESPremixed2D_pState &d_dWdy_dy) const;

   /**************** Operators Overloading ********************/
   /* Index operator */
   double &operator[](int index);
   const double &operator[](int index) const;

   /* Binary arithmetic operators. */
   LESPremixed2D_pState operator +(const LESPremixed2D_pState &W) const;
   LESPremixed2D_pState operator -(const LESPremixed2D_pState &W) const;
   LESPremixed2D_pState operator *(const double &a) const;
   friend LESPremixed2D_pState operator *(const double &a, const LESPremixed2D_pState &W);
   LESPremixed2D_pState operator /(const double &a) const;
   double operator *(const LESPremixed2D_pState &W) const;
   LESPremixed2D_pState operator ^(const LESPremixed2D_pState &W) const;

   /* Assignment Operator. */ 
   LESPremixed2D_pState& operator =(const LESPremixed2D_pState &W); 

   /* Shortcut arithmetic operators. */
   LESPremixed2D_pState& operator +=(const LESPremixed2D_pState &W);
   LESPremixed2D_pState& operator -=(const LESPremixed2D_pState &W);
 
   /* Unary arithmetic operators. */
   //LESPremixed2D_pState operator +(const LESPremixed2D_pState &W);
   friend LESPremixed2D_pState operator -(const LESPremixed2D_pState &W);

   /* Relational operators. */
   friend int operator ==(const LESPremixed2D_pState &W1, const LESPremixed2D_pState &W2);
   friend int operator !=(const LESPremixed2D_pState &W1, const LESPremixed2D_pState &W2);

   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, const LESPremixed2D_pState &W);
   friend istream& operator >> (istream &in_file,  LESPremixed2D_pState &W);

   /**************** Destructors ******************************/
   //for static specdata
   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; 
                                           specdata = NULL; 
				 if(Schmidt != NULL) delete[] Schmidt; 
 				           Schmidt = NULL;
                                }

   //for scalars
   void Deallocate_addvar(void){ if(scalar != NULL)  delete[] scalar;   
                                           scalar = NULL;
                               }

   void scalar_memory() { Deallocate_addvar(); 
                          if (nscal) scalar = new double[nscal];
                        }

#ifdef STATIC_LESPREMIXED2D_SPECIES    
   void specnull() {}
   void spec_memory() {}
   void Deallocate_species(void) {}
   void Deallocate(void){ Deallocate_addvar(); }
#else
   void specnull() {spec=NULL;}
   void spec_memory() { Deallocate_species(); spec = new Species[ns];}
   void Deallocate_species(void) { if(spec != NULL){ delete[] spec;}  specnull(); }
   void Deallocate(void){ Deallocate_species(); Deallocate_addvar(); }
#endif

   ~LESPremixed2D_pState(){ Deallocate(); }

};


 /**************************************************************************
 ********************* LESPREMIXED2D_CSTATE CLASS DEFINTION ****************
     The conserved variables edition of LESPREMIXED2D_PSTATE... 
 ***************************************************************************
 ***************************************************************************/
 class LESPremixed2D_cState {
   private: 
   //all public .... yes I know ....
   protected:
   public:
   /****** Conservative Vector "U" **************************/
  double                      rho; //!< Density.
  Vector2D                   rhov; //!< Momentum (2D)
  double                        E; //!< Total Energy (rho *(e + HALF*v^2))
#ifdef STATIC_LESPREMIXED2D_SPECIES
   Species   rhospec[STATIC_LESPREMIXED2D_SPECIES];
#else 
  Species                *rhospec;    //!< Species class using (rho*c[ns])
#endif
  double                *rhoscalar;   //!< Scalars using (rho*scalar[nscal])

   Tensor2D                    tau;   //!< Shear Stress Tensor
   Vector2D                  qflux;   //!< Heat Flux Vector  
   Tensor2D                 lambda;   //!< SFS Stress Tensor
   Vector2D                  theta;   //!< Turbulent Heat Flux Vector  

#ifdef THICKENED_FLAME_ON
   PowerLaw                  flame; //!< SFS wrinkling factor and thickening factor
#endif 
 
   static int                         ns; //!< number of species
   static int                      nscal; //!< number of scalars
   static NASARP1311data       *specdata; //!< Global species data 
   static double                *Schmidt; //!< Schmidt Number for each species
   static Set_scalar            Scal_sys; //!< Set the group of scalars to be solved in the model
   static double          low_temp_range; //!< Low temp data range
   static double         high_temp_range; //!< High temp data range
   static int      NUM_VAR_LESPREMIXED2D; //!< Number of LESPremixed2D variables (4+ns)
   static double                    Mref; //!< Mref for Precondtioning (normally set to incoming freestream Mach)
   static SubfilterScaleModels  SFSmodel; //!< Subfilter scale model
   static double            filter_width; //!< Constant filter width
   //@}
 
   //@{ @name Premixed flame parameters:
   static double      laminar_speed; //!< Unstrained laminar flame speed 
   static double  laminar_thickness; //!< Unstrained laminar flame thickness
   static double            TFactor; //!< Maximum thickening factor
   static double     adiabatic_temp; //!< Adiabatic temperature
   static double  equivalence_ratio; //!< Equivalence ratio
   static double      reactants_den; //!< Reactants density
   //@}



   /*****************************************************/
   //default constructors of many flavours, hopefully one is right 4U
   LESPremixed2D_cState(): rho(DENSITY_STDATM), E(PRESSURE_STDATM/(rho*((double)0.4)))
                    { rhoscalar = NULL; set_initial_values_scal();  rhospecnull(); set_initial_values(); }   

   LESPremixed2D_cState(const double &value): rho(value), rhov(value), E(value)
                 { rhoscalar = NULL; set_initial_values_scal();  rhospecnull();  set_initial_values(value); }   

   LESPremixed2D_cState(const double &d, const Vector2D &V, const double &En):
                 rho(d), rhov(V), E(En)
 		{ rhoscalar = NULL; set_initial_values_scal(); rhospecnull(); set_initial_values(); }

   LESPremixed2D_cState(const double &d, const double &vx, const double &vy, const double &En):
                 rho(d), rhov(vx,vy), E(En)
 		{ rhoscalar = NULL; set_initial_values_scal();  rhospecnull(); set_initial_values(); }

   LESPremixed2D_cState(const double &d, const Vector2D &V, const double &En, const double *rhoscal):
                 rho(d), rhov(V), E(En)
 		{ rhoscalar = NULL; set_initial_values_scal(rhoscal);  rhospecnull(); set_initial_values(); }

   LESPremixed2D_cState(const double &d, const Vector2D &V, const double &En, const double &value):
                 rho(d), rhov(V), E(En)
		 { rhoscalar = NULL; set_initial_values_scal();  rhospecnull(); set_initial_values(value); }

   LESPremixed2D_cState(const double &d, const double &vx, const double &vy, const double &En, 
		 const double *rhoscal):
                 rho(d), rhov(vx,vy), E(En)
     { rhoscalar = NULL; set_initial_values_scal(rhoscal);  rhospecnull();   set_initial_values(); }

   LESPremixed2D_cState(const double &d, const double &vx, const double &vy, const double &En, const double &value): 
                 rho(d), rhov(vx,vy), E(En)
                 { rhoscalar = NULL; set_initial_values_scal();  rhospecnull();  set_initial_values(value); }

   LESPremixed2D_cState(const double &d, const double &vx, const double &vy, const double &En, 
                 const double *rhoscal, const double &value): 
                 rho(d), rhov(vx,vy), E(En)
 		{ rhoscalar = NULL; set_initial_values_scal(rhoscal);  rhospecnull();   set_initial_values(value); }

   LESPremixed2D_cState(const double &d, const double &vx, const double &vy, const double &En,
 		const double *rhoscal, const Species *rhomfrac): 
                 rho(d), rhov(vx,vy), E(En)
 		{ rhoscalar = NULL; set_initial_values_scal(rhoscal); rhospecnull();   set_initial_values(rhomfrac); }

   LESPremixed2D_cState(const double d, const Vector2D &V, const double &En, const double *rhoscal, const Species *rhomfrac): 
                 rho(d), rhov(V), E(En)
     { rhoscalar = NULL; set_initial_values_scal(rhoscal);  rhospecnull();   set_initial_values(rhomfrac); }

   // WARNING - automatic type conversion
   LESPremixed2D_cState(const LESPremixed2D_pState &W) : rho(W.rho), rhov(W.rhov()), E(W.E()), 
		   tau(W.tau), qflux(W.qflux), lambda(W.lambda), theta(W.theta) 
#ifdef THICKENED_FLAME_ON
		   ,flame(W.flame)
#endif		   
   {
     for(int i=0; i<W.ns; ++i){
       rhospec[i].c = W.rho*W.spec[i].c;
       rhospec[i].gradc = W.rho*W.spec[i].gradc;
       rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
     }
     if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] = W.rho*W.scalar[i];    
   }

   //this is needed for the operator overload returns!!!!
   LESPremixed2D_cState(const LESPremixed2D_cState &U): rho(U.rho), rhov(U.rhov), E(U.E),
		                  tau(U.tau), qflux(U.qflux), lambda(U.lambda), theta(U.theta)
#ifdef THICKENED_FLAME_ON
		   ,flame(U.flame)
#endif
		                  { rhoscalar = NULL; set_initial_values_scal(U.rhoscalar); 
				    rhospecnull(); set_initial_values(U.rhospec); }


   //read in ns species data, call only once as its static
   void set_species_data(const int &, const int &, const string *, const char *,
			 const double&, const double *, const int &);

   //set initial data values predominately used internally 
   void set_initial_values();
   void set_initial_values(const double &value);
   void set_initial_values(double *rhomfrac);
   void set_initial_values(const Species *rhomfrac);

   void set_initial_values_scal();
   void set_initial_values_scal(const double &value);
   void set_initial_values_scal(const double *rhoscal);


   //Copy construtor, cheaper than = operator
   void Copy(const LESPremixed2D_cState &U);

   /***************** VACUUM ************************/
   void Vacuum(){ rho=ZERO; rhov.zero(); E=ZERO;  
     if (nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] = ZERO;
     for(int i=0; i<ns; ++i) rhospec[i].Vacuum();
     tau.zero();  qflux.zero(); lambda.zero(); theta.zero(); 
   }  

   void zero_non_sol(){
     for(int i=0; i<ns; ++i){
       rhospec[i].gradc.zero();
       rhospec[i].diffusion_coef=ZERO;
     }
     tau.zero(); qflux.zero(); lambda.zero(); theta.zero();  
   }  


   //! Set premixed flame static variables.
   void set_premixed_flame_variables(const double &lam_thickness,
				     const double &lam_speed,
				     const double &thickening_factor,
                                     const double &adia_temp,
                                     const double &equi_ratio,
                                     const double &reac_den);
   //! Set SFS modelling variables.
   void set_SFSmodel_variables(const double &delta,
			       const double &Smagorinsky_coef,
			       const double &Yoshizawa_coef);
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
   double k(void) const;
   double T(void) const;      //temperature
   double a(void) const;      //speed of sound
   double rhok(void) const; 
   /************  FSD Model*******************/
   double C(void) const;
   double Fsd(void) const;
   bool negative_scalarcheck(void) const;

   bool negative_speccheck(const int &step) ; //-ve mass frac check and sets small -ve c's to ZERO
   bool Unphysical_Properties();
   bool Unphysical_Properties_Check(const int Flow_Type, const int n);   
   double sum_species(void) const;

   //@{ @name Functions required for multigrid.
   //! Copy variables solved by multigrid only.
   void Copy_Multigrid_State_Variables(const LESPremixed2D_cState &Ufine);
  
   //! Zero variables not-solved by multigrid.
   void Zero_Non_Multigrid_State_Variables(void);
   //@}


   /************** Temperature Derivatives *******************/
   double dmudT(void) const;

   // STRAIN RATE, LAMINAR STRESS, SFS STRESS 6 FUNCTIONS\
   // LIKE PRIMITIVE 
   /************ Strain rate tensor, laminar stress and SFS stress tensors ***********/
   Tensor2D Strain_Rate(const LESPremixed2D_pState &dWdx,
			const LESPremixed2D_pState &dWdy,
			const int Flow_Type,
			const int Axisymmetric,
			const Vector2D X) const;

   /************ Heat Flux vector thermal Diffusion ***********/
   Vector2D thermal_diffusion(const double &Temp) const;

   /*********** Primitive solution state ***********************/
   LESPremixed2D_pState W(void) const;
   LESPremixed2D_pState W(const LESPremixed2D_cState &U) const;
   friend LESPremixed2D_pState W(const LESPremixed2D_cState &U);

   /**************** turbulence model related parameters*********/
   double mu_t(const Tensor2D &strain_rate, const int &Flow_Type) const;      
   double Pr_turb(void) const;      
   double Sc_turb(void) const;      
   double Dm_turb(const double &mut) const;  
  

   /**************** Fluxes ***********************************/
   //Viscous Flux (laminar+turbulent)
   LESPremixed2D_cState Viscous_Flux_x(const LESPremixed2D_pState &dWdx, 
				       const LESPremixed2D_pState &dWdy, 
				       const int Flow_Type,
				       const int Axisymmetric,
				       const Vector2D X) const;
   LESPremixed2D_cState Viscous_Flux_y(const LESPremixed2D_pState &dWdx, 
				       const LESPremixed2D_pState &dWdy, 
				       const int Flow_Type,
				       const int Axisymmetric,
				       const Vector2D X) const;

   /*************** Preconditioner ****************************/
   void Low_Mach_Number_Preconditioner(DenseMatrix &P,const int &flow_type_flag, const double &deltax) const; 
   void Low_Mach_Number_Preconditioner_Inverse(DenseMatrix &Pinv,const int &flow_type_flag, const double &deltax) const; 

   /************** Premixed combustion ************************/
   LESPremixed2D_cState  premixed_mfrac(const LESPremixed2D_pState &Wo);
   int FlameJump_x(const LESPremixed2D_cState &Uu,
		   LESPremixed2D_cState &Ub,
		   const int &Flow_Type);
   int FlameJump_n(const LESPremixed2D_cState &Uu,
		   LESPremixed2D_cState &Ub,
		   const Vector2D &norm_dir,
		   const int &Flow_Type);

   /***************** Index operators *************************/
   double &operator[](int index);
   const double &operator[](int index) const;

   /**************** Operators Overloading ********************/
   /* Binary arithmetic operators. */
   LESPremixed2D_cState operator +(const LESPremixed2D_cState &U) const;
   LESPremixed2D_cState operator -(const LESPremixed2D_cState &U) const;
   LESPremixed2D_cState operator *(const double &a) const;
   friend LESPremixed2D_cState operator *(const double &a, const LESPremixed2D_cState &U);
   LESPremixed2D_cState operator /(const double &a) const;

   double operator *(const LESPremixed2D_cState &U) const;
   LESPremixed2D_cState operator ^(const LESPremixed2D_cState &U) const;

   /* Assignment Operator. */ 
   LESPremixed2D_cState& operator =(const LESPremixed2D_cState &U); 

   /* Shortcut arithmetic operators. */
   LESPremixed2D_cState& operator +=(const LESPremixed2D_cState &U);
   LESPremixed2D_cState& operator -=(const LESPremixed2D_cState &U);
   LESPremixed2D_cState& operator *=(const double &a);
   LESPremixed2D_cState& operator /=(const double &a);

   /* Unary arithmetic operators. */
   //LESPremixed2D_cState operator +(const LESPremixed2D_cState &U);
   friend LESPremixed2D_cState operator -(const LESPremixed2D_cState &U);

   /* Relational operators. */
   friend int operator ==(const LESPremixed2D_cState &U1, const LESPremixed2D_cState &U2);
   friend int operator !=(const LESPremixed2D_cState &U1, const LESPremixed2D_cState &U2);

   /* Input-output operators. */
   friend ostream& operator << (ostream &out_file, const LESPremixed2D_cState &U);
   friend istream& operator >> (istream &in_file,  LESPremixed2D_cState &U);

   /**************** Destructors ******************************/
   void Deallocate_static(void){ if(specdata != NULL) delete[] specdata; 
                                   specdata = NULL; 
 				 if(Schmidt != NULL) delete[] Schmidt; 
 				   Schmidt = NULL; 
                               }

   //for scalars and additional variables
   void Deallocate_addvar(void){ if(rhoscalar != NULL) delete [] rhoscalar;   
                                   rhoscalar = NULL; 
                               }

   void rhoscalar_memory() { Deallocate_addvar(); 
                             if (nscal) rhoscalar = new double[nscal];
                           }

#ifdef STATIC_LESPREMIXED2D_SPECIES  
   void rhospecnull() {}          
   void rhospec_memory() {}
   void Deallocate_species() {}
   void Deallocate(void) { Deallocate_addvar(); }
#else
   void rhospecnull() {rhospec=NULL;}
   void rhospec_memory() { Deallocate_species(); rhospec = new Species[ns];}
   void Deallocate_species() { if(rhospec != NULL) { delete[] rhospec; } rhospecnull();}
   void Deallocate(void){ Deallocate_species(); Deallocate_addvar(); }   
#endif

   ~LESPremixed2D_cState() { Deallocate(); }

 };


 /**************************************************************************
 ****************** LESPREMIXED2D_PSTATE CONSTRUCTORS **********************
 ***************************************************************************/

 /**************************************************************************
    Set up mass fractions memory and initial values, ie. Species Class  
 **************************************************************************/
 inline void LESPremixed2D_pState::set_initial_values(){
   spec_memory();
   for(int i=0; i<ns; ++i) spec[i].c = ONE/ns; 
 }

 inline void  LESPremixed2D_pState::set_initial_values(const double &value){
   spec_memory();
   for(int i=0; i<ns; ++i) spec[i].c = value; 
 }

 //user specified
 inline void LESPremixed2D_pState::set_initial_values(double *cfrac){
   spec_memory();
   for(int i=0; i<ns; ++i) spec[i].c = cfrac[i];
 }

 //another set using species class
 inline void LESPremixed2D_pState::set_initial_values(const Species *mfrac){
   spec_memory();
   for(int i=0; i<ns; ++i) spec[i] = mfrac[i];
 }

/**************************************************************************
   Set up scalars memory and initial values  
**************************************************************************/
inline void LESPremixed2D_pState::set_initial_values_scal(){
  scalar_memory();
  for(int i=0; i<nscal; ++i) scalar[i] = ZERO; 
}

inline void LESPremixed2D_pState::set_initial_values_scal(const double &value){
  scalar_memory();
  for(int i=0; i<nscal; ++i) scalar[i] = value; 
}

//user specified
inline void LESPremixed2D_pState::set_initial_values_scal(const double *scal){
  scalar_memory();
  for(int i=0; i<nscal; ++i) scalar[i] = scal[i];  
}

/***************************************************************************
 * LESPremixed2D_pState::set_SFSmodel_variables -- Set the SFS modelling   *
 *                                                 variables.              *
 ***************************************************************************/
inline void LESPremixed2D_pState::set_SFSmodel_variables(const double &delta,
							 const double &Smagorinsky_coef,
							 const double &Yoshizawa_coef) {
  filter_width = delta;
  SFSmodel.set_SFSmodel_constants(Smagorinsky_coef, Yoshizawa_coef);
}

/***************************************************************************
 * LESPremixed2D_pState::set_premixed_flame_variables -- Set the premixed  *
 *                                                flame static variables.  *
 ***************************************************************************/
inline void LESPremixed2D_pState::set_premixed_flame_variables(const double &lam_thickness,
							       const double &lam_speed,
							       const double &thickening_factor,
                                                               const double &adia_temp,
                                                               const double &equi_ratio,
                                                               const double &reac_den){
  laminar_speed = lam_speed; 
  laminar_thickness = lam_thickness;
  TFactor = thickening_factor; 
  adiabatic_temp = adia_temp;
  equivalence_ratio = equi_ratio;
  reactants_den = reac_den;
}


 /*****************  Momentum *******************************/
 inline Vector2D LESPremixed2D_pState::rhov(void) const{
   return rho*v; 
 } 

 /********************** Prandtl ****************************/
 inline double LESPremixed2D_pState::Prandtl(void) const{
   //Pr = Cp*mu/k
   return Cp()*mu()/kappa();
 }

 /********************** Schmidt ****************************/
 inline double LESPremixed2D_pState::Schmidt_No(const int &i) const{
   if(spec[i].diffusion_coef > ZERO){
#ifdef THICKENED_FLAME_ON
     return (mu()/(flame.WF*flame.TF))/(rho*spec[i].diffusion_coef);
#else
     return mu()/(rho*spec[i].diffusion_coef);
#endif
   } else {
     return Schmidt[i];
   }
 }

 /********************** Lewis *****************************/
 inline double LESPremixed2D_pState::Lewis(const int &i) const{
   if(spec[i].diffusion_coef > ZERO){
#ifdef THICKENED_FLAME_ON
     return (kappa()/(flame.WF*flame.TF))/(rho*Cp()*spec[i].diffusion_coef);
#else
     return kappa()/(rho*Cp()*spec[i].diffusion_coef);
#endif
   }
   return ZERO;
 }

 /******* Mixture Diffusion Coefficient ********************/
 // inline double LESPremixed2D_pState::Diffusion_coef(void) const{
 //   double sum=ZERO;
 //   for(int i=0; i<ns; ++i){
 //     sum += spec[i].c * spec[i].diffusion_coef;
 //   }
 //   return sum;
 // }

 /************* Temperature ********************************/
 inline double LESPremixed2D_pState::T(void) const{
   return p/(rho*Rtot());
 }

 /************* Turbulence modified  pressure **************/
 inline double LESPremixed2D_pState::pmodified(void) const{
   return ( p + 2.0*rho*k()/3.0 );
 }

 /************** SFS turbulence kinetic energy *************/
 inline double LESPremixed2D_pState::k(void) const{
   return (-lambda.trace())/(TWO*rho);
 }

 /************* Strain rate tensor ***********************************/
 inline Tensor2D LESPremixed2D_pState::Strain_Rate(const LESPremixed2D_pState &dWdx,
						   const LESPremixed2D_pState &dWdy,
						   const int Flow_Type,
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
 inline void LESPremixed2D_pState::Laminar_Stress(const Tensor2D &strain_rate){
   tau = TWO*mu()*strain_rate;
 }

 /***************** Turbulent SFS stress ******************************/
inline void LESPremixed2D_pState::SFS_Stress(const Tensor2D &strain_rate,
                                             const int &Flow_Type){

  double kk = SFSmodel.sfs_k_Yoshizawa(strain_rate, filter_width);

   lambda = TWO * mu_t(strain_rate,Flow_Type) * strain_rate;
   lambda.xx -= (TWO/THREE)*rho*kk;
   lambda.yy -= (TWO/THREE)*rho*kk;
   lambda.zz -= (TWO/THREE)*rho*kk;
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
inline bool LESPremixed2D_pState::negative_speccheck(void) {
  double sum(ZERO);
  //double SPEC_TOLERANCE(SPEC_TOLERANCE);

  //-------- Negative Check ------------//
  for(int i=0; i<ns-1; ++i){
    if(spec[i].c > ONE){ //check for > 1.0
      spec[i].c = ONE;
    } else if(spec[i].c < ZERO){  //check for -ve
      if(spec[i].c > -SPEC_TOLERANCE){  //check for small -ve and set to ZERO 
	spec[i].c = ZERO;
      } else {
	spec[i].c = ZERO;
	//#ifdef _DEBUG
	cout <<"\n pState -ve mass fraction in "<<specdata[i].Speciesname()<<" "<<
	  spec[i].c<<" greater than allowed tolerance of "<<-SPEC_TOLERANCE; 
	//#endif
      }      
    } 
    sum += spec[i].c;
  } 
  
  // spec[ns-1].c = (ONE - sum);  //PUSH error into NS-1
  // Spread error across species
  spec[ns-1].c = max(ONE- sum, ZERO);
  sum += spec[ns-1].c;
  for(int i=0; i<ns; ++i){
    spec[i].c = spec[i].c*(ONE/sum);
  }

  // check for scalars
  if (nscal) {
    for (int i=0; i<nscal; ++i) if(scalar[i] < 0.0) scalar[i] = 0.0; 
  }

  return true;
}

/***** Species Concentrations ******************************/
inline double LESPremixed2D_pState::SpecCon(int i) const{
  //returned in kg/m^3 / kg/mol => mol/m^3
  return (rho)*spec[i].c/(specdata[i].Mol_mass());
}

/******* GIBBS Free Energy ********************************
  Eqn. (10.84) (10.87) Anderson
  Gs = Hs - TS
  Gs(ps=1) = Gs - R_UNIVERSAL*T*ln(ps)  //for data not at 1atm
  ps = cs(M/Ms)*p
***********************************************************/
inline double LESPremixed2D_pState::Gibbs(int species) const{
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
inline double LESPremixed2D_pState::diedip() const{
  double Rlocal = Rtot();
  return (hprime() - Rlocal)/(rho*Rlocal);
}

inline double LESPremixed2D_pState::diedirho() const{
  double Rlocal = Rtot();
  return -p*(hprime() - Rlocal)/(rho*rho*Rlocal);
}

/**************** Copy *************************************/
inline void LESPremixed2D_pState::Copy(const LESPremixed2D_pState &W){
  rho = W.rho;
  v = W.v; 
  p = W.p;  
  if(nscal) for(int i =0; i<nscal; ++i) scalar[i] = W.scalar[i];
  for( int i=0; i<ns; ++i) spec[i] = W.spec[i];
  tau = W.tau;
  qflux = W.qflux;
  lambda = W.lambda;
  theta = W.theta;
#ifdef THICKENED_FLAME_ON
  flame = W.flame;
#endif

}

//**************** Index Operators *************************/
inline double& LESPremixed2D_pState::operator[](int index) {  
  //  assert( index >= 1 && index <= NUM_VAR_LESPREMIXED2D );

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
  if(index <= NUM_VAR_LESPREMIXED2D-ns) {
      return scalar[index-NUM_LESPREMIXED2D_VAR_SANS_SPECIES-1];
  } else if(index <= NUM_VAR_LESPREMIXED2D) {
      return spec[index-(NUM_VAR_LESPREMIXED2D-ns)-1].c;
  }
#ifdef THICKENED_FLAME_ON
      if(index == NUM_VAR_LESPREMIXED2D+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
#endif
   } 
};

inline const double& LESPremixed2D_pState::operator[](int index) const {  
  //  assert( index >= 1 && index <= NUM_VAR_LESPREMIXED2D );

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
  if(index <= NUM_VAR_LESPREMIXED2D-ns) {
      return scalar[index-NUM_LESPREMIXED2D_VAR_SANS_SPECIES-1];
  } else if(index <= NUM_VAR_LESPREMIXED2D) {
      return spec[index-(NUM_VAR_LESPREMIXED2D-ns)-1].c;
  }
#ifdef THICKENED_FLAME_ON
      if (index == NUM_VAR_LESPREMIXED2D+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
#endif
    }      
};

/**************************************************************
  Get max of the min temperature of the lowest region
  and min of the max temperature of the highest region
***************************************************************/
inline void LESPremixed2D_pState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; ++i){
    temp = max(specdata[i].Low_range(),temp);
  }
  low_temp_range = temp;  
}

inline void LESPremixed2D_pState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; ++i){
    temp = min(specdata[i].High_range(),temp);
  }
  high_temp_range = temp;  
}

/********************************************************
 * LESPremixed2D_pState::U -- Conserved solution state. *
 ********************************************************/
inline LESPremixed2D_cState LESPremixed2D_pState::U(void) const {
  return U(*this);
}

inline LESPremixed2D_cState LESPremixed2D_pState::U(const LESPremixed2D_pState &W) const{
    LESPremixed2D_cState Temp;
    Temp.rho = W.rho;
    Temp.rhov = W.rhov();
    for(int i=0; i<W.ns; ++i){
      Temp.rhospec[i].c = W.rho*W.spec[i].c;
      Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
      Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
    }
    if(nscal) for(int i=0; i<nscal; ++i) Temp.rhoscalar[i] = W.rho*W.scalar[i];    
    Temp.E = W.E();
    
    Temp.tau = W.tau;
    Temp.qflux = W.qflux; 
    Temp.lambda = W.lambda;
    Temp.theta = W.theta; 
#ifdef THICKENED_FLAME_ON
    Temp.flame = W.flame;
#endif

    return Temp;
}

inline LESPremixed2D_cState U(const LESPremixed2D_pState &W) {
  LESPremixed2D_cState Temp;
  Temp.rho = W.rho;
  Temp.rhov = W.rhov();
  for(int i=0; i<W.ns; ++i){
    Temp.rhospec[i].c = W.rho*W.spec[i].c;
    Temp.rhospec[i].gradc = W.rho*W.spec[i].gradc;
    Temp.rhospec[i].diffusion_coef = W.rho*W.spec[i].diffusion_coef;
  }  
  if(W.nscal) for(int i=0; i<W.nscal; ++i) Temp.rhoscalar[i] = W.rho*W.scalar[i];    
  Temp.E = W.E();
  
  Temp.tau = W.tau;
  Temp.qflux = W.qflux; 
  Temp.lambda = W.lambda;
  Temp.theta = W.theta; 
#ifdef THICKENED_FLAME_ON
  Temp.flame = W.flame;
#endif

  return Temp;
}

/**************************************************************************
***************** LESPREMIXED2D_CSTATE CONSTRUCTORS ***********************
***************************************************************************/

/***************************************************************
  Set up mass fractions memory and initial values, ie. Species Class  
****************************************************************/
inline void LESPremixed2D_cState::set_initial_values(){
  rhospec_memory();
  for(int i=0; i<ns; ++i) rhospec[i].c = rho/ns; 
} 

inline void  LESPremixed2D_cState::set_initial_values(const double &value){
  rhospec_memory();
  for(int i=0; i<ns; ++i) rhospec[i].c = value; 
}

//user specified
inline void LESPremixed2D_cState::set_initial_values(double *rhomfrac){
  rhospec_memory();
  for(int i=0; i<ns; ++i) rhospec[i].c = rhomfrac[i];
}

//another set using species class
inline void LESPremixed2D_cState::set_initial_values( const Species *rhomfrac){
  rhospec_memory();
  for(int i=0; i<ns; ++i) rhospec[i] = rhomfrac[i];
}

/**************************************************************************
   Set up scalars memory and initial values
**************************************************************************/
inline void LESPremixed2D_cState::set_initial_values_scal(){
  rhoscalar_memory();
  for(int i=0; i<nscal; ++i) rhoscalar[i] = ZERO; 
}

inline void LESPremixed2D_cState::set_initial_values_scal(const double &value){
  rhoscalar_memory();
  for(int i=0; i<nscal; ++i) rhoscalar[i] = value;  
}

//user specified
inline void LESPremixed2D_cState::set_initial_values_scal(const double *rhoscal){
  rhoscalar_memory();
  for(int i=0; i<nscal; ++i) rhoscalar[i] = rhoscal[i];     
}

/**************** Copy *************************************/
inline void LESPremixed2D_cState::Copy(const LESPremixed2D_cState &U){
  rho = U.rho;
  rhov = U.rhov; 
  E = U.E;
  if(nscal) for(int i=0; i<nscal; ++i) rhoscalar[i] = U.rhoscalar[i];
  for( int i=0; i<ns; ++i) rhospec[i] = U.rhospec[i]; 
  tau = U.tau;
  qflux = U.qflux; 
  lambda = U.lambda;
  theta = U.theta;
#ifdef THICKENED_FLAME_ON
  flame = U.flame;
#endif    
}

/***************************************************************************
 * LESPremixed2D_cState::set_SFSmodel_variables -- Set the SFS modelling   *
 *                                                 variables.              *
 ***************************************************************************/
inline void LESPremixed2D_cState::set_SFSmodel_variables(const double &delta,
							 const double &Smagorinsky_coef,
							 const double &Yoshizawa_coef) {
  filter_width = delta;
  SFSmodel.set_SFSmodel_constants(Smagorinsky_coef, Yoshizawa_coef);
}

/***************************************************************************
 * LESPremixed2D_cState::set_premixed_flame_variables -- Set the premixed  *
 *                                                flame static variables.  *
 ***************************************************************************/
inline void LESPremixed2D_cState::set_premixed_flame_variables(const double &lam_thickness,
							       const double &lam_speed,
							       const double &thickening_factor,
                                                               const double &adia_temp,
                                                               const double &equi_ratio,
                                                               const double &reac_den){
  laminar_speed = lam_speed; 
  laminar_thickness = lam_thickness;
  TFactor = thickening_factor; 
  adiabatic_temp = adia_temp;
  equivalence_ratio = equi_ratio;
  reactants_den = reac_den;
}


/**************** Velocity *********************************/
inline Vector2D LESPremixed2D_cState::v() const{
  return (rhov/rho);
}

/**************** Pressure *********************************/
inline double LESPremixed2D_cState::p() const{
  return (rho*Rtot()*T());
}
/**************** kinetic energy ***************************/
inline double LESPremixed2D_cState::k() const{
  return rhok()/rho;
}

inline double LESPremixed2D_cState::rhok() const{
  return (-lambda.trace())/TWO;
}

/***************** C ************************/
inline double LESPremixed2D_cState::C(void) const{
  if(nscal > 0){
    return rhoscalar[0]/rho;
  } else {
    return 0.0;
  }
}

/**************** Fsd ************************/
inline double LESPremixed2D_cState::Fsd(void) const{
  if(nscal > 0){
    return rhoscalar[1]/rho;
  } else { 
    return 0.0;
  }
}

/************* Strain rate tensor ***********************************/
inline Tensor2D LESPremixed2D_cState::Strain_Rate(const LESPremixed2D_pState &dWdx,
						  const LESPremixed2D_pState &dWdy,
						  const int Flow_Type,
						  const int Axisymmetric,
						  const Vector2D X) const{
   Tensor2D strain_rate;
   double radius, div_v;

   /***************** Strain rate (+dilatation) **********************/	
   div_v = dWdx.v.x + dWdy.v.y;
   if (Axisymmetric == AXISYMMETRIC_X) {
     radius = (X.x < MICRO) ? MICRO : X.x;    //fabs(X.x) ??
     div_v += v().x/radius;
   } else if (Axisymmetric == AXISYMMETRIC_Y) {    
     radius = (X.y < MICRO) ? MICRO : X.y;
     div_v += v().y/radius;
   } 

   strain_rate.xx = dWdx.v.x-div_v/THREE;
   strain_rate.xy = HALF*(dWdx.v.y + dWdy.v.x);
   strain_rate.yy = dWdy.v.y-div_v/THREE;

   if (Axisymmetric == PLANAR) {
     strain_rate.zz = -(strain_rate.xx + strain_rate.yy); 
   } else if (Axisymmetric == AXISYMMETRIC_X) {
     strain_rate.zz = v().x/radius-div_v/THREE;
   } else if (Axisymmetric == AXISYMMETRIC_Y) {
     strain_rate.zz = v().y/radius-div_v/THREE;
   } 

   return strain_rate;  
}


//----------------- Index Operator ---------------------//
inline double& LESPremixed2D_cState::operator[](int index) {
  //  assert( index >= 1 && index <= NUM_VAR_LESPREMIXED2D );
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
  if(index <= NUM_VAR_LESPREMIXED2D-ns) {
      return rhoscalar[index-NUM_LESPREMIXED2D_VAR_SANS_SPECIES-1];
  } else if(index <= NUM_VAR_LESPREMIXED2D) {
      return rhospec[index-(NUM_VAR_LESPREMIXED2D-ns)-1].c;
  }
#ifdef THICKENED_FLAME_ON
      if (index == NUM_VAR_LESPREMIXED2D+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
#endif
    }
};

inline const double& LESPremixed2D_cState::operator[](int index) const{
  //  assert( index >= 1 && index <= NUM_VAR_LESPREMIXED2D ); 
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
   if(index <= NUM_VAR_LESPREMIXED2D-ns) {
      return rhoscalar[index-NUM_LESPREMIXED2D_VAR_SANS_SPECIES-1];
  } else if(index <= NUM_VAR_LESPREMIXED2D) {
      return rhospec[index-(NUM_VAR_LESPREMIXED2D-ns)-1].c;
  }
#ifdef THICKENED_FLAME_ON
      if (index == NUM_VAR_LESPREMIXED2D+1) {
	return flame.WF;
      } else {
	return flame.TF;
      }
#endif
    }
};

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
inline bool LESPremixed2D_cState::negative_speccheck(const int &step) {
  double sum(ZERO);
  double temp(ZERO);

  //-------- Negative Check ------------//     
  for(int i=0; i<ns-1; ++i){
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
  for(int i=0; i<ns; ++i){
    rhospec[i].c = rhospec[i].c*(ONE/sum);
  } 


  // check for scalars
  if (nscal) {
    for (int i=0; i<nscal; ++i) if(rhoscalar[i] < 0.0) rhoscalar[i] = 0.0; 
  }
  
  return true;
}

// USED IN MULTIGRID
inline bool LESPremixed2D_cState::Unphysical_Properties(){
  LESPremixed2D_cState::Unphysical_Properties_Check(FLOWTYPE_LAMINAR,10);
}

inline bool LESPremixed2D_cState::negative_scalarcheck(void) const{
  double LOCAL_TOL = MICRO; 
  //-------- Negative Check ------------//
  for ( int i=0; i<nscal; i++ ) {
    if ( rhoscalar[i] < -LOCAL_TOL ) {
      rhoscalar[i] = ZERO;
    }
  }
      if ( nscal > 1 ) {
        if ( rhoscalar[0]/rho > 0.99 || rhoscalar[0]/rho < 0.01 ) {
        rhoscalar[1] = ZERO;
      }
      }
  return (1);
}

/**************************************************************
  Unphysical_Properties_Check
 ***************************************************************/
inline bool LESPremixed2D_cState::Unphysical_Properties_Check(const int Flow_Type, const int n){

  // check for nan's, inf's etc.... debugging !!!
#ifdef _DEBUG
  for( int i=0; i<NUM_LESPREMIXED2D_VAR_SANS_SPECIES+ns; ++i){ 
    if( this[i] != this[i]){ cout<<"\n nan's in solution, variable "<<i<<endl; exit(1); return false; }
  }
#endif

  if ((Flow_Type == FLOWTYPE_INVISCID ||
       Flow_Type == FLOWTYPE_LAMINAR) &&
      (rho <= ZERO || !negative_speccheck(n) || es() <= ZERO)) {
    cout << "\n " << CFFC_Name() 
	 << " LESPremixed2D ERROR: Negative Density || Energy || mass fractions: \n" << *this <<endl;
    return false;
  }
  if ((Flow_Type == FLOWTYPE_TURBULENT_LES_TF_SMAGORINSKY || 
       Flow_Type == FLOWTYPE_TURBULENT_LES_TF_K) &&
      (rho <= ZERO || !negative_speccheck(n) || es() <= ZERO || rhok() < ZERO)) {
    cout << "\n " << CFFC_Name() 
	 << " LESPremixed2D ERROR: Negative Density || Energy || mass fractions || Turbulent kinetic energy || : \n"
	 << *this <<endl;
    return false;
  }
 if ((Flow_Type == FLOWTYPE_LAMINAR_C ||                            
      Flow_Type == FLOWTYPE_LAMINAR_C_ALGEBRAIC ||       
      Flow_Type == FLOWTYPE_LAMINAR_C_FSD || 
      Flow_Type == FLOWTYPE_LAMINAR_NGT_C_FSD || 
      Flow_Type == FLOWTYPE_TURBULENT_LES_C ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_C_ALGEBRAIC || 
      Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_SMAGORINSKY ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_CHARLETTE ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_NGT_C_FSD_SMAGORINSKY ||
      Flow_Type == FLOWTYPE_TURBULENT_LES_C_FSD_K ||
      Flow_Type == FLOWTYPE_FROZEN_TURBULENT_LES_C_FSD) &&
     (rho <= ZERO || !negative_speccheck(n) || es() <= ZERO || !negative_scalarcheck() )) {
    cout << "\n " << CFFC_Name() 
	 << " LESPremixed2D ERROR: Negative Density || Energy || Mass Fractions || Progress Variable || Flame Surface Density : \n"
	 << *this << endl;
    return false;
  }else{
    return true ;
  }
} 

/**************************************************************
  Sum N-1 species.
***************************************************************/
inline double LESPremixed2D_cState::sum_species(void) const{
  double sum=ZERO;
  for(int i=0; i<ns-1; ++i){
    sum += rhospec[i].c;
  }
  return sum/rho;
}

/**********************************************************************
 * LESPremixed2D_cState::Copy_Multigrid_State_Variables --            *
 *                           Copy variables solved by multigrid only. *
 **********************************************************************/
inline void LESPremixed2D_cState::Copy_Multigrid_State_Variables(const LESPremixed2D_cState &Ufine) {
  Copy(Ufine);
  Zero_Non_Multigrid_State_Variables();
} 
  
/**********************************************************************
 * LESPremixed2D_cState::Zero_Non_Multigrid_State_Variables --        *
 *                            Zero variables not-solved by multigrid. *
 **********************************************************************/
inline void LESPremixed2D_cState::Zero_Non_Multigrid_State_Variables(void) {
  for(int i=0; i<nscal; ++i) rhoscalar[i] = ZERO;
}

/**************************************************************
  Get max and min temperature ranges for data
***************************************************************/
inline void LESPremixed2D_cState::Temp_low_range(void){  
  double temp = specdata[0].Low_range();
  for(int i=0; i<ns; ++i){
    temp = max(specdata[i].Low_range(),temp);
  } 
  low_temp_range = temp;  
}

inline void LESPremixed2D_cState::Temp_high_range(void){
  double temp = specdata[0].High_range();
  for(int i=0; i<ns; ++i){
    temp = min(specdata[i].High_range(),temp);
  } 
  high_temp_range = temp;  
}


/********************************************************
 * LESPremixed2D_cState::W -- Primitive solution state. *
 ********************************************************/
inline LESPremixed2D_pState LESPremixed2D_cState::W(void) const {
  return W(*this);
}

inline LESPremixed2D_pState LESPremixed2D_cState::W(const LESPremixed2D_cState &U) const{
    LESPremixed2D_pState Temp;
    Temp.rho = U.rho;
    Temp.v = U.v();  
    for(int i=0; i<U.ns; ++i){
      Temp.spec[i].c = U.rhospec[i].c/U.rho;
      Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
      Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
    }   
    Temp.p = U.p();
    if(nscal) for(int i=0; i<nscal; ++i) Temp.scalar[i] = U.rhoscalar[i]/U.rho;
    
    Temp.tau = U.tau;
    Temp.qflux = U.qflux; 
    Temp.lambda = U.lambda;
    Temp.theta = U.theta; 
#ifdef THICKENED_FLAME_ON
    Temp.flame = U.flame;
#endif   

    return Temp;
}

inline LESPremixed2D_pState W(const LESPremixed2D_cState &U) {
  LESPremixed2D_pState Temp;
  Temp.rho = U.rho;
  Temp.v = U.v();
  for(int i=0; i<U.ns; ++i){
    Temp.spec[i].c = U.rhospec[i].c/U.rho;
    Temp.spec[i].gradc = U.rhospec[i].gradc/U.rho;
    Temp.spec[i].diffusion_coef = U.rhospec[i].diffusion_coef/U.rho;
  } 
  Temp.p = U.p();
  if(U.nscal) for(int i=0; i<U.nscal; ++i) Temp.scalar[i] = U.rhoscalar[i]/U.rho;
    
  Temp.tau = U.tau;
  Temp.qflux = U.qflux;
  Temp.lambda = U.lambda;
  Temp.theta = U.theta;
#ifdef THICKENED_FLAME_ON
  Temp.flame = U.flame;
#endif   

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
extern LESPremixed2D_pState Reflect(const LESPremixed2D_pState &W,
				    const Vector2D &norm_dir);

extern LESPremixed2D_pState Free_Slip(const LESPremixed2D_pState &Win,
				      const LESPremixed2D_pState &Wout,
				      const Vector2D &norm_dir,
				      const int &TEMPERATURE_BC_FLAG);

extern LESPremixed2D_pState No_Slip(const LESPremixed2D_pState &Win,
				    const LESPremixed2D_pState &Wout,
				    const Vector2D &norm_dir,
				    const int &TEMPERATURE_BC_FLAG);

extern LESPremixed2D_pState Moving_Wall(const LESPremixed2D_pState &Win,
					const LESPremixed2D_pState &Wout,
					const Vector2D &norm_dir,				 
					const double &wall_velocity,
					const int &TEMPERATURE_BC_FLAG);

extern LESPremixed2D_pState BC_Characteristic_Pressure(const LESPremixed2D_pState &Wi,
						       const LESPremixed2D_pState &Wo,
						       const Vector2D &norm_dir);

extern LESPremixed2D_pState BC_1DFlame_Inflow(const LESPremixed2D_pState &Wi,
					      const LESPremixed2D_pState &Wo, 
					      const LESPremixed2D_pState &Woutlet,
					      const Vector2D &norm_dir);

extern LESPremixed2D_pState BC_2DFlame_Inflow(const LESPremixed2D_pState &Wi,
					      const LESPremixed2D_pState &Wo, 				
					      const Vector2D &norm_dir);

extern LESPremixed2D_pState BC_1DFlame_Outflow(const LESPremixed2D_pState &Wi,
					       const LESPremixed2D_pState &Wo,
					       const LESPremixed2D_pState &Winlet,
					       const Vector2D &norm_dir);

extern LESPremixed2D_pState BC_2DFlame_Outflow(const LESPremixed2D_pState &Wi, 
					       const LESPremixed2D_pState &Wo,
					       const Vector2D &norm_dir);

/*******************************************************
 * Exact Test Case Solution Functions                  *
 *******************************************************/
extern LESPremixed2D_pState RinglebFlow(const LESPremixed2D_pState &Wdum,
					const Vector2D X);

extern LESPremixed2D_pState ViscousChannelFlow(const LESPremixed2D_pState &Wdum,
					       const Vector2D X,
					       const double Vwall,
					       const double dp);

extern LESPremixed2D_pState FlatPlate(const LESPremixed2D_pState &Winf,
				      const Vector2D X,
				      double &eta,
				      double &f,
				      double &fp,
				      double &fpp);

/*******************************************************
 * External Flux Function Functions                    *
 *******************************************************/
extern LESPremixed2D_pState RoeAverage(const LESPremixed2D_pState &Wl,
				       const LESPremixed2D_pState &Wr);

// HLLE
extern LESPremixed2D_cState FluxHLLE_x(const LESPremixed2D_pState &Wl,
				       const LESPremixed2D_pState &Wr,
				       const int &Preconditioning, 
				       const int Flow_Type);

extern LESPremixed2D_cState FluxHLLE_x(const LESPremixed2D_cState &Ul,
				       const LESPremixed2D_cState &Ur,
				       const int &Preconditioning,
				       const int Flow_Type);
  
extern LESPremixed2D_cState FluxHLLE_n(const LESPremixed2D_pState &Wl,
				       const LESPremixed2D_pState &Wr,
				       const Vector2D &norm_dir,
				       const int &Preconditioning,
				       const int Flow_Type);

extern LESPremixed2D_cState FluxHLLE_n(const LESPremixed2D_cState &Ul,
				       const LESPremixed2D_cState &Ur,
				       const Vector2D &norm_dir,
				       const int &Preconditioning);

// Linde
extern LESPremixed2D_cState FluxLinde(const LESPremixed2D_pState &Wl,
				      const LESPremixed2D_pState &Wr,
				      const int Flow_Type);

extern LESPremixed2D_cState FluxLinde(const LESPremixed2D_cState &Ul,
				      const LESPremixed2D_cState &Ur,
				      const int Flow_Type);


extern LESPremixed2D_cState FluxLinde_n(const LESPremixed2D_pState &Wl,
					const LESPremixed2D_pState &Wr,
					const Vector2D &norm_dir,
					const int Flow_Type);

extern LESPremixed2D_cState FluxLinde_n(const LESPremixed2D_cState &Ul,
					const LESPremixed2D_cState &Ur,
					const Vector2D &norm_dir,
					const int Flow_Type);

// Roe
extern LESPremixed2D_pState WaveSpeedPos(const LESPremixed2D_pState &lambda_a,
					 const LESPremixed2D_pState &lambda_l,
					 const LESPremixed2D_pState &lambda_r);
			

extern LESPremixed2D_pState WaveSpeedNeg(const LESPremixed2D_pState &lambda_a,
					 const LESPremixed2D_pState &lambda_l,
					 const LESPremixed2D_pState &lambda_r);

extern LESPremixed2D_pState WaveSpeedAbs(const LESPremixed2D_pState &lambda_a,
					 const LESPremixed2D_pState &lambda_l,
					 const LESPremixed2D_pState &lambda_r);
 
extern LESPremixed2D_pState HartenFixPos(const LESPremixed2D_pState &lambda_a,
					 const LESPremixed2D_pState &lambda_l,
					 const LESPremixed2D_pState &lambda_r );
				

extern LESPremixed2D_pState HartenFixNeg(const LESPremixed2D_pState &lambda_a,
					 const LESPremixed2D_pState &lambda_l,
					 const LESPremixed2D_pState &lambda_r);
				  

extern LESPremixed2D_pState HartenFixAbs(const LESPremixed2D_pState &lambdas_a,
					 const LESPremixed2D_pState &lambdas_l,
					 const LESPremixed2D_pState &lambdas_r );
				

extern LESPremixed2D_cState FluxRoe_x(const LESPremixed2D_pState &Wl,
				      const LESPremixed2D_pState &Wr,
				      const int &Preconditioning, 
				      const int &flow_type_flag,
				      const double &deltax);

extern LESPremixed2D_cState FluxRoe_x(const LESPremixed2D_cState &Ul,
				      const LESPremixed2D_cState &Ur,
				      const int &Preconditioning, 
				      const int &flow_type_flag,
				      const double &deltax);

extern LESPremixed2D_cState FluxRoe_n(const LESPremixed2D_pState &Wl,
				      const LESPremixed2D_pState &Wr,
				      const Vector2D &norm_dir,
				      const int &Preconditioning, 
				      const int &flow_type_flag,
				      const double &deltax);

extern LESPremixed2D_cState FluxRoe_n(const LESPremixed2D_cState &Ul,
				      const LESPremixed2D_cState &Ur,
				      const Vector2D &norm_dir,
				      const int &Preconditioning, 
				      const int &flow_type_flag,
				      const double &deltax);

extern LESPremixed2D_cState FluxAUSMplus_up(const LESPremixed2D_pState &Wl,
					    const LESPremixed2D_pState &Wr);

extern LESPremixed2D_cState FluxAUSMplus_up(const LESPremixed2D_cState &Ul,
					    const LESPremixed2D_cState &Ur);

extern LESPremixed2D_cState FluxAUSMplus_up_n(const LESPremixed2D_pState &Wl,
					      const LESPremixed2D_pState &Wr,
					      const Vector2D &norm_dir);

extern LESPremixed2D_cState FluxAUSMplus_up_n(const LESPremixed2D_cState &Ul,
					      const LESPremixed2D_cState &Ur,
					      const Vector2D &norm_dir);

/* Viscous Solution flux (laminar+turbulent) */
extern LESPremixed2D_cState Viscous_FluxArithmetic_n(const LESPremixed2D_cState &Ul,
						     const LESPremixed2D_pState &dWdx_l,
						     const LESPremixed2D_pState &dWdy_l,
						     const LESPremixed2D_cState &Ur,
						     const LESPremixed2D_pState &dWdx_r,
						     const LESPremixed2D_pState &dWdy_r,
						     const int Flow_Type,
						     const int Axisymmetric,
						     const Vector2D X,
						     const Vector2D &norm_dir);

extern LESPremixed2D_cState Viscous_Flux_n(LESPremixed2D_pState &W,
					   const LESPremixed2D_pState &dWdx,
					   const LESPremixed2D_pState &dWdy,
					   const int Flow_Type,
					   const int Axisymmetric,
					   const Vector2D X,			
					   const Vector2D &norm_dir);

extern double WallShearStress(const LESPremixed2D_pState &W1,
			      const Vector2D &X1,
			      const Vector2D &X2,
			      const Vector2D &X3,
			      const Vector2D &norm_dir);

/*********************************************************/
extern LESPremixed2D_pState Rotate(const LESPremixed2D_pState &W,
				   const Vector2D &norm_dir);


extern Vector2D HLLE_wavespeeds(const LESPremixed2D_pState &Wl,
                                const LESPremixed2D_pState &Wr,
                                const Vector2D &norm_dir); 

#endif //end _LESPREMIXED2D_STATE_INCLUDED 
