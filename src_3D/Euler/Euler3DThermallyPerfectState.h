/*! \file Euler3DThermallyPerfectState.h
 * 	\brief	Header file defining the Euler solution state classes 
 *              associated with solution of compressible inviscid flows 
 *              of a thermally perfect non-reactive or combusting mixture.
 *
 *  \note  If #define STATIC_NUMBER_OF_SPECIES is set dynamic memory
 *         is not used so code is faster, however code is not as flexible
 *         as only up to STATIC_NUMBER_OF_SPECIES can be used without 
 *         recompling the entire code.
 */

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

/* Define the classes. */

class Euler3D_ThermallyPerfect_cState;
class Euler3D_ThermallyPerfect_pState;

/*! Number of fixed variables in the Euler3D_ThermallyPerfect class */
#define NUM_EULER3D_VAR_SANS_SPECIES 5  //rho, v(3), p

/*! If you define this variable, the number of species will be 
    predetermined for faster calculations.., however it is not as general */
#define STATIC_NUMBER_OF_SPECIES 5 // set to 2 for air, 5 for 1-step CH4, 6 for 2-setp CH4

/*! Convergence tolerance */
#define CONV_TOLERANCE 1e-8

/*!
 * Class: Euler3D_ThermallyPerfect_pState
 *
 * \brief Primitive state solution class for 3D Euler equations
 *        governing flows of thermally perfect non-reactive and 
 *        combusting mixtures.
 *
 * Member functions
 *  - rho             -- Return density (kg/m^3)
 *  - v               -- Return flow velocity (m/s)
 *  - p               -- Return pressure (Pa, N/m^2)
 *  - spec            -- Return array of species mass fraction data
 *  - Mass            -- Return mixture molecular mass (kg/mol)
 *  - Rtot            -- Return mixture gas constant (J/(kg*K))
 *  - HeatofFormation -- Return heat of formation for the mixture
 *  - Cp              -- Return specific heat at constant pressure for mixture (J/(kg*K))
 *  - Cv              -- Return specific heat at constant volume for mixture (J/(kg*K))
 *  - g               -- Return specific heat ratio for mixture
 *  - e               -- Return mixture absolute internal energy (J/kg)
 *  - es              -- Return mixture sensible internal energy (J/kg)
 *  - h               -- Return mixture absolute specific enthalpy (J/kg)
 *  - hs              -- Return mixture sensible specific enthalpy (J/kg)
 *  - E               -- Return total mixture energy (J/kg)
 *  - H               -- Return total mixture enthalpy (J/kg)
 *  - rhov            -- Return momentum of mixture (kg/(m^2*s))
 *  - T               -- Return mixture temperature (K)
 *  - a               -- Return sound speed of mixture (m/s)
 *  - M               -- Return Mach number for mixture
 *  - Gibbs           -- Return species Gibbs free energy
 *  - U               -- Return conserved solution state
 *  - F, Fx           -- Return x-direction inviscid solution flux
 *  - Fy              -- Return y-direction inviscid solution flux
 *  - Fz              -- Return z-direction inviscid solution flux
 *  - Schem           -- Return source terms associated with finite-rate chemistry
 *  - lambda          -- Return x-direction eigenvalue
 *  - rc              -- Return x-direction conserved right eigenvector
 *  - lp              -- Return x-direction primitive left eigenvector
 *  - lambda_x        -- Return x-direction eigenvalue
 *  - rc_x            -- Return x-direction conserved right eigenvector
 *  - lp_x            -- Return x-direction primitive left eigenvector
 *  - RoeAverage      -- Return Roe-average solution state vector
 *  - FluxHLLE_x      -- Return HLLE numerical solution flux in x-direction
 *  - FluxHLLE_n      -- Return HLLE numerical solution flux in n-direction
 *  - FluxRoe_x       -- Return Roe numerical solution flux in x-direction
 *  - FluxRoe_n       -- Return Roe numerical solution flux in n-direction
 *  - lambda_minus    -- Return negative eigenvalues, applying Harten entropy fix
 *  - lambda_plus     -- Return positive eigenvalues, applying Harten entropy fix
 *  - Reflect         -- Return reflected solution state after application of reflection BC
 *  - MovingWall      -- Return wall solution state after application of moving wall BC
 *  - NoSlip          -- Return wall solution state after application of no-slip BC
 *
 * Member operators \n
 *      W -- a primitive solution state \n
 *      c -- a scalar (double)
 *
 *  - W = W;
 *  - c = W[i];
 *  - W = W + W;
 *  - W = W - W;
 *  - c = W * W; (inner product)
 *  - W = c * W;
 *  - W = W * c;
 *  - W = W / c;
 *  - W = W ^ W;
 *  - W = +W;
 *  - W = -W;
 *  - W += W;
 *  - W -= W;
 *  - W *= W;
 *  - W /= W;
 *  - W == W;
 *  - W != W;
 *  - cout << W; (output function)
 *  - cin  >> W; (input function)
 *  \nosubgrouping
 */
class Euler3D_ThermallyPerfect_pState {
  protected:
   static double             *_temp_values;  //!< Temporary static values for species calculations
   static double              *_diff_coeff;  //!< Temporary static values for species diffusion calculations

  public:
   double                              rho;  //!< Mixture density (kg/m^3)
   Vector3D                              v;  //!< Mixture velocity vector (m/s)
   double                                p;  //!< Mixture pressure (Pa)

#ifdef STATIC_NUMBER_OF_SPECIES
   Species  spec[STATIC_NUMBER_OF_SPECIES];  //!< Species concentration data for multi-component mixture
#else 
   Species                           *spec;  //!< Species concentration data for multi-component mixture 
#endif

   static int                     num_vars;  //!< Number of solution variables or unknowns
   static int                           ns;  //!< Number of species in multi-component mixture
   static NASARP1311data         *specdata;  //!< Species transport data for multi-component mixture   
   static double                  *Schmidt;  //!< Species Schmidt numbers
   static Reaction_set               React;  //!< Reaction mechanism data
   static double            low_temp_range;  //!< Lower temperature bound for thermodynamic data
   static double           high_temp_range;  //!< Upper temperature bound for thermodynamic data   
   static int                  debug_level;  //!< Level for debug information
   
/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   Euler3D_ThermallyPerfect_pState(void): 
     rho(DENSITY_STDATM), p(PRESSURE_STDATM) {
     v.zero(); species_null(); set_initial_values();
   }
   
   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &value) : 
     rho(value), p(value) {
     v.x = value;  v.y = value; v.z = value; 
     species_null(); set_initial_values(value);
   }
         
   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V,
                                   const double &pre) :
     rho(d), v(V), p(pre) {
     species_null(); set_initial_values();
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V, 
                                   const double &pre, 
                                   const double &frac):
     rho(d), v(V), p(pre) {
     species_null(); set_initial_values(frac);
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V, 
                                   const double &pre, 
                                   const Species *mfrac) :
     rho(d), v(V), p(pre) {
     species_null(); set_initial_values(mfrac); 
   }
  
   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const Vector3D &V,
	 	                   const double &pre, 
                                   const double *mfrac) :
     rho(d), v(V), p(pre) {
     species_null(); set_initial_values(mfrac); 
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx,
                                   const double &vy, 
                                   const double &vz,
                                   const double &pre) :
     rho(d), p(pre) {
     v.x = vx;  v.y = vy; v.z = vz;
     species_null(); set_initial_values(); 
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &pre, 
                                   const double &frac) :
     rho(d), p(pre) {
     v.x = vx; v.y = vy; v.z = vz; 
     species_null(); set_initial_values(frac); 
   }
                    
   //! Assignment constructor
   Euler3D_ThermallyPerfect_pState(const double &d, 
                                   const double &vx, 
                                   const double &vy,  
                                   const double &vz,
                                   const double &pre, 
                                   const Species *mfrac) :
     rho(d), p (pre) {
     v.x=vx; v.y=vy; v.z = vz; 
     species_null(); set_initial_values(mfrac); 
   }
  
   //! Copy constructor (this is needed for the operator overload returns)
   Euler3D_ThermallyPerfect_pState(const Euler3D_ThermallyPerfect_pState &W) {
     species_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(W);
   }
     
   //! Default destructor
   ~Euler3D_ThermallyPerfect_pState(void) {
      Deallocate(); 
   }
//@}

/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Allocate memory for species data
   void Allocate(void) {
#ifdef STATIC_NUMBER_OF_SPECIES
      if (STATIC_NUMBER_OF_SPECIES < ns) {
        cerr << "\n WARNING USING STATIC ARRAYS FOR SPECIES DATA WITH STATIC_NUMBER_OF_SPECIES = " 
             << STATIC_NUMBER_OF_SPECIES << "." 
             << "\n HOWEVER THERE IS A REQUEST FOR FOR ns = " << ns << " SPECIES.\n"; 
        exit(1); 
      }
#else
      if (spec != NULL) {
         delete[] spec;
         spec = NULL;
      } /* endif */
      spec = new Species[ns];
#endif
   }

   //! Allocate memory for species data
   void Allocate(const int &n) {
#ifdef STATIC_NUMBER_OF_SPECIES
      if (STATIC_NUMBER_OF_SPECIES < n) {
        cerr << "\n WARNING USING STATIC ARRAYS FOR SPECIES DATA WITH STATIC_NUMBER_OF_SPECIES = " 
             << STATIC_NUMBER_OF_SPECIES << "." 
             << "\n HOWEVER THERE IS A REQUEST FOR FOR n = " << n << " SPECIES.\n"; 
        exit(1); 
      }
#else
      if (spec != NULL) {
         delete[] spec;
         spec = NULL;
      } /* endif */
      ns = n;
      spec = new Species[ns];
#endif
   }

   //! Deallocate memory for species data
   void Deallocate(void) {
#ifndef STATIC_NUMBER_OF_SPECIES 
      if (spec != NULL) delete[] spec;
      spec = NULL;
#endif
   }

   //! Allocate static memory for species data
   void Allocate_static(void) {
      if (_temp_values != NULL) {
         delete[] _temp_values;
         _temp_values = NULL;
      } /* endif */
      _temp_values = new double[ns];
      if (_diff_coeff != NULL) {
         delete[] _diff_coeff;
         _diff_coeff = NULL;
      } /* endif */
      _diff_coeff = new double[ns];
      if (specdata != NULL) {
         delete[] specdata;
         specdata = NULL;
      } /* endif */
      specdata = new NASARP1311data[ns]; 
      if (Schmidt != NULL) {
         delete[] Schmidt; 
         Schmidt = NULL;
      } /* endif */
      Schmidt = new double[ns];
   }

   //! Allocate static memory for species data
   void Allocate_static(const int &n) {
      ns = n;
      if (_temp_values != NULL) {
         delete[] _temp_values;
         _temp_values = NULL;
      } /* endif */
      _temp_values = new double[ns];
      if (_diff_coeff != NULL) {
         delete[] _diff_coeff;
         _diff_coeff = NULL;
      } /* endif */
      _diff_coeff = new double[ns];
      if (specdata != NULL) {
         delete[] specdata;
         specdata = NULL;
      } /* endif */
      specdata = new NASARP1311data[ns]; 
      if (Schmidt != NULL) {
         delete[] Schmidt; 
         Schmidt = NULL;
      } /* endif */
      Schmidt = new double[ns];
   }

   //! Deallocate static memory for species data
   void Deallocate_static(void) {
      if (_temp_values != NULL) delete[] _temp_values;
      _temp_values = NULL;
      if (_diff_coeff != NULL) delete[] _diff_coeff;
      _diff_coeff = NULL;
      if (specdata != NULL) delete[] specdata;
      specdata = NULL;
      if (Schmidt != NULL) delete[] Schmidt; 
      Schmidt = NULL;
   }

   //! Set species point to null
   void species_null(void) {
#ifndef STATIC_NUMBER_OF_SPECIES 
      spec = NULL;
#endif
   }

   //! Assign the species data (needs to be called only once as it is static)
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);
     
   //! Set initial values for species data
   void set_initial_values(void);
   //! Set initial values for species data
   void set_initial_values(const double &value);
   //! Set initial values for species data
   void set_initial_values(const double *cfrac);
   //! Set initial values for species data
   void set_initial_values(const Species *mfrac);

   //! Copy solution state (cheaper than = operator)
   void Copy(const Euler3D_ThermallyPerfect_pState &W);

   //! Returns the number of variables - number of species
   int NumVarSansSpecies(void) const { return num_vars - ns; }

   //! Assigns a vacuum solution state
   void Vacuum(void){ 
     rho=ZERO; v.x=ZERO; v.y=ZERO; v.z=ZERO;  p = ZERO;
     for (int i = 0; i < ns; i++) {
        spec[i].Vacuum();
     } /* endfor */
   }
  
   //! Set the lower bound for valid temperature ranges for thermodynamic data 
   void Temp_low_range(void);     
   //! Set the upper bound for valid temperature ranges for thermodynamic data 
   void Temp_high_range(void); 

  //! Check for negative species mass fractions
   bool negative_speccheck(void);

   //! Check for physical validity of the solution vector
   bool Realizable_Solution_Check(void);
 //@}

/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
   //! Mixture molecular mass
   double Mass(void) const;
   
   //! Mixture gas constant
   double Rtot(void);        
   //! Mixture gas constant
   double Rtot(void) const;   

   //! Mixture heat of formation
   double HeatofFormation(void);        
   //! Mixture heat of formation
   double HeatofFormation(void) const;

   //! Mixture specific heat at constant pressure
   double Cp(void) const;    
   //! Mixture specific heat at constant pressure
   double Cp(const double& TEMP) const;

   //! Mixture specific heat at constant volume
   double Cv(void) const;     

   //! Mixture specific heat ratio
   double g(void) const;    
   //! Estimate for mixture specific heat ratio
   double gamma_guess(void) const; 

   //! Mixture absolute (sensible+chemical) internal energy
   double e(void) const;     

   //! Mixture sensible internal energy  
   double es(void) const; 

   //! Reference internal energy (sensible + heat of formation - offset)
   double eref(void) const;

   //! Mixture absolute specific enthalpy   
   double h(void) const;      
   //! Mixture absolute specific enthalpy
   double h(const double &T) const;

   //! Mixture sensible specific enthalpy
   double hs(void) const;
   //! Mixture sensible specific enthalpy
   double hs(const double &T) const;

   //! Reference internal enthalpy (sensible + heat of formation - offset)
   double href(void) const;

   //! Total energy of mixture
   double E(void) const;      

   //! Total enthalpy of mixture
   double H(void) const;

   //! Total sensible enthalpy of mixture
   double Hs(void) const;    

   //! Derivative of specific enthalpy wrt temperature, dh/dT
   double hprime(void) const;  
   //! Derivative of specific enthalpy wrt temperature, dh/dT
   double hprime(double &Temp) const;

   //! Mixture momentum
   Vector3D rhov(void) const;      

   //! Mixture temperature
   double T(void) const;         
   //! Mixture temperature given sensible enthalpy 
   double T(double &h_s) const;   

   //! Mixture sound speed
   double a(void);                
   //! Mixture sound speed
   double a(void) const;

   //! Mixture Mach number
   double M(void) const;
     
   //! Gibbs Free Energy (H-TS) for species
   double Gibbs(int species) const;

   //! Species i concentration (rho*c/mol_mass)
   double SpecCon(int i) const;

   //! Partial derivative of energy wrt pressure
   double diedip(void) const;

   //! Partial derivative of energy wrt density
   double diedirho(void) const;
//@}

/** @name Conserved solution state */ 
/*        ------------------------ */
//@{
   //! Returns conserved solution state
   Euler3D_ThermallyPerfect_cState U(void); 
   //! Returns conserved solution state
   Euler3D_ThermallyPerfect_cState U(void) const; 
   //! Returns conserved solution state
   Euler3D_ThermallyPerfect_cState U(const Euler3D_ThermallyPerfect_pState &W); 
//@}
 
/** @name Inviscid Flux Vectors */
/*        --------------------- */
//@{
   //! x-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState F(void);
   //! x-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState F(void) const;

   //! x-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fx(void);
   //! x-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fx(void) const;

   //! y-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fy(void);
   //! y-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fy(void) const;

   //! z-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fz(void);
   //! z-direction inviscid solution flux
   Euler3D_ThermallyPerfect_cState Fz(void) const;
//@}

/** @name Flux Jacobians */
/*        -------------- */
//@{
   //!< x-direction inviscid flux Jacobian
   void dFxdU(DenseMatrix &dFxdU);
//@}

/** @name Eigenvalue(s) and eigenvector(s) (x-direction) */
/*        ---------------------------------------------- */
//@{
   //! x-direction eigenvalues
   Euler3D_ThermallyPerfect_pState lambda(void);
   //! x-direction eigenvalues
   Euler3D_ThermallyPerfect_pState lambda(void) const;

   //! x-direction eigenvalues
   Euler3D_ThermallyPerfect_pState lambda_x(void);
   //! x-direction eigenvalues
   Euler3D_ThermallyPerfect_pState lambda_x(void) const;

   //! x-direction conservative eigenvectors
   Euler3D_ThermallyPerfect_cState rc(const int &index); 
   //! x-direction conservative eigenvectors
   Euler3D_ThermallyPerfect_cState rc(const int &index) const;

   //! x-direction conservative eigenvectors
   Euler3D_ThermallyPerfect_cState rc_x(const int &index); 
   //! x-direction conservative eigenvectors
   Euler3D_ThermallyPerfect_cState rc_x(const int &index) const;

   //! x-direction primitive eignenvectors
   Euler3D_ThermallyPerfect_pState lp(const int &index) ; 
   //! x-direction primitive eignenvectors
   Euler3D_ThermallyPerfect_pState lp(const int &index) const;

   //! x-direction primitive eignenvectors
   Euler3D_ThermallyPerfect_pState lp_x(const int &index) ; 
   //! x-direction primitive eignenvectors
   Euler3D_ThermallyPerfect_pState lp_x(const int &index) const; 
//@}
  
/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns Roe-averaged primitive solution state
   static Euler3D_ThermallyPerfect_pState RoeAverage(const Euler3D_ThermallyPerfect_pState &Wl,
                                                     const Euler3D_ThermallyPerfect_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static Euler3D_ThermallyPerfect_cState FluxHLLE_x(const Euler3D_ThermallyPerfect_pState &Wl,
                                                     const Euler3D_ThermallyPerfect_pState &Wr);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) x-direction flux
   static Euler3D_ThermallyPerfect_cState FluxHLLE_x(const Euler3D_ThermallyPerfect_cState &Ul,
                                                     const Euler3D_ThermallyPerfect_cState &Ur);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static Euler3D_ThermallyPerfect_cState FluxHLLE_n(const Euler3D_ThermallyPerfect_pState &Wl,
                                                     const Euler3D_ThermallyPerfect_pState &Wr,
                                                     const Vector3D &norm_dir);

   //! Returns Harten-Lax-van-Leer-Einfeldt (HLLE) flux in n-direction
   static Euler3D_ThermallyPerfect_cState FluxHLLE_n(const Euler3D_ThermallyPerfect_cState &Ul,
                                                     const Euler3D_ThermallyPerfect_cState &Ur,
                                                     const Vector3D &norm_dir);

   //! Returns Roe flux x-direction flux
   static Euler3D_ThermallyPerfect_cState FluxRoe_x(const Euler3D_ThermallyPerfect_pState &Wl, 
                                                    const Euler3D_ThermallyPerfect_pState &Wr);

   //! Returns Roe flux in n-direction
   static Euler3D_ThermallyPerfect_cState FluxRoe_n(const Euler3D_ThermallyPerfect_pState &Wl,
                                                    const Euler3D_ThermallyPerfect_pState &Wr,
                                                    const Vector3D &norm_dir);

   //! Returns negative waves speeds (eigenvalues) using Harten entropy fix
   static Euler3D_ThermallyPerfect_pState lambda_minus(const Euler3D_ThermallyPerfect_pState  &lambda_a,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_l,
                                                       const Euler3D_ThermallyPerfect_pState  &lambda_r);

   //! Returns positive waves speeds (eigenvalues) using Harten entropy fix
   static Euler3D_ThermallyPerfect_pState lambda_plus(const Euler3D_ThermallyPerfect_pState  &lambda_a,
                                                      const Euler3D_ThermallyPerfect_pState  &lambda_l,
                                                      const Euler3D_ThermallyPerfect_pState  &lambda_r); 

   //! HLLE wavespeeds in n-direction given 2 primitive states and a direction
   static Vector2D HLLE_wavespeeds(const Euler3D_ThermallyPerfect_pState &Wl,
                                   const Euler3D_ThermallyPerfect_pState &Wr,
                                   const Vector3D &norm_dir);
  
   //! Returns rotated primitive state aligned with local x-axis in the norm_dir
   Euler3D_ThermallyPerfect_pState Rotate(const Vector3D &norm_dir) const;				  
  
   //! Returns un-rotated primitive state aligned with x-axis of global problem
   Euler3D_ThermallyPerfect_pState RotateBack(const Vector3D &norm_dir) const;
//@}
				  
/** @name Boundary Conditions */
/*        ------------------- */
//@{
   //! Return reflected solution state after application of reflection BC
   static Euler3D_ThermallyPerfect_pState Reflect(const Euler3D_ThermallyPerfect_pState &W, 
                                                  const Vector3D &norm_dir);

   //! Return wall solution state after application of moving wall BC
   static Euler3D_ThermallyPerfect_pState MovingWall(const Euler3D_ThermallyPerfect_pState &Win,
                                                     const Euler3D_ThermallyPerfect_pState &Wout,
                                                     const Vector3D &norm_dir,
                                                     const Vector3D &wall_velocity,
                                                     const Vector3D &pressure_gradient,
                                                     const int &TEMPERATURE_BC_FLAG);

   //! Return wall solution state after application of no-slip BC
   static Euler3D_ThermallyPerfect_pState NoSlip(const Euler3D_ThermallyPerfect_pState &Win, 
                                                 const Euler3D_ThermallyPerfect_pState &Wout,
                                                 const Vector3D &norm_dir,  
                                                 const Vector3D &pressure_gradient,
                                                 const int &TEMPERATURE_BC_FLAG);
//@}

/** @name Finite-rate chemical kinetics */
/*        ----------------------------- */
//@{
   //! Returns source solution vector associated with finite-rate chemical kinetics
   Euler3D_ThermallyPerfect_cState Schemistry(int &REACT_SET_FLAG) const;

   //! Returns Jacobian of chemical reaction rate source vector wrt conserved solution vector
   void dSchemistrydU(DenseMatrix &dSwdU) const;

   //! Returns maximum absolute value of diagonal entry of Jacobian of source vector
   double dSchemistrydU_max_diagonal(void) const;
//@}  

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;

   //! Binary addition operator
   Euler3D_ThermallyPerfect_pState operator +(const Euler3D_ThermallyPerfect_pState &W) const;

   //! Binary subtraction operator
   Euler3D_ThermallyPerfect_pState operator -(const Euler3D_ThermallyPerfect_pState &W) const;

   //! Binary multiplication operator
   Euler3D_ThermallyPerfect_pState operator *(const double &a) const;

   //! Binary multiplication operator
   friend Euler3D_ThermallyPerfect_pState operator *(const double &a, 
                                                     const Euler3D_ThermallyPerfect_pState &W);

   //! Binary multiplication operator
   double operator *(const Euler3D_ThermallyPerfect_pState &W) const;

   //! Binary division operator
   Euler3D_ThermallyPerfect_pState operator /(const double &a) const;

   //! Binary vector product operator
   Euler3D_ThermallyPerfect_pState operator ^(const Euler3D_ThermallyPerfect_pState &W) const;
  
   //! Assignment operator
   Euler3D_ThermallyPerfect_pState& operator =(const Euler3D_ThermallyPerfect_pState &W); 
  
   //! Shortcut addition operator
   Euler3D_ThermallyPerfect_pState& operator +=(const Euler3D_ThermallyPerfect_pState &W);

   //! Shortcut subtraction operator
   Euler3D_ThermallyPerfect_pState& operator -=(const Euler3D_ThermallyPerfect_pState &W);

   //! Unary subtraction operators
   friend Euler3D_ThermallyPerfect_pState operator -(const Euler3D_ThermallyPerfect_pState &W);
  
   //! Equal relational operator
   friend int operator ==(const Euler3D_ThermallyPerfect_pState &W1, 
                          const Euler3D_ThermallyPerfect_pState &W2);

   //! Not equal relational operator
   friend int operator !=(const Euler3D_ThermallyPerfect_pState &W1, 
                          const Euler3D_ThermallyPerfect_pState &W2);

   //! Output stream operator
   friend ostream& operator << (ostream &out_file, 
                                const Euler3D_ThermallyPerfect_pState &W);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,  
                                Euler3D_ThermallyPerfect_pState &W);
//@}
};

/*!
 * Class: Euler3D_ThermallyPerfect_cState
 *
 * \brief Conserved state solution class for 3D Euler equations
 *        governing flows of thermally perfect non-reactive and 
 *        combusting mixtures.
 *
 * Member functions
 *  - rho     -- Return mixture density (kg/m^3)
 *  - rhov    -- Return mixture momentum (kg/(m^2-s))   
 *  - E       -- Return mixture total energy (J/kg)
 *  - rhospec -- Return array of species density data
 *
 * Member operators \n
 *      U -- a conserved solution state \n
 *      c -- a scalar (double)
 *
 *  - U = U;
 *  - c = U[i];
 *  - U = U + U;
 *  - U = U - U;
 *  - c = U * U; (inner product)
 *  - U = c * U;
 *  - U = U * c;
 *  - U = U / c;
 *  - U = U ^ U;
 *  - U = +U;
 *  - U = -U;
 *  - U += U;
 *  - U -= U;
 *  - U *= U;
 *  - U /= U;
 *  - U == U;
 *  - U != U;
 *  - cout << U; (output function)
 *  - cin  >> U; (input function)
 *  \nosubgrouping
 */
class Euler3D_ThermallyPerfect_cState {
  protected:
   static double                *_temp_values;  //!< Temporary static values for species calculations
   static double                 *_diff_coeff;  //!< Temporary static values for species diffusion calculations

  public:
   double                                 rho;  //!< Mixture density (kg/m^3) 
   Vector3D                              rhov;  //!< Mixture momentum (kg/(m^2-s))    
   double                                   E;  //!< Mixture total energy (J/kg)

#ifdef STATIC_NUMBER_OF_SPECIES
   Species  rhospec[STATIC_NUMBER_OF_SPECIES];  //!< Species density data for multi-component mixture
#else 
   Species                           *rhospec;  //!< Species density data for multi-component mixture
#endif

   static int                        num_vars;  //!< Number of solution variables or unknowns 
   static int                              ns;  //!< Number of species in multi-component mixture
   static NASARP1311data            *specdata;  //!< Species transport data for multi-component mixture
   static double                     *Schmidt;  //!< Species Schmidt numbers
   static double               low_temp_range;  //!< Lower temperature bound for thermodynamic data
   static double              high_temp_range;  //!< Upper temperature bound for thermodynamic data 
   static int                     debug_level;  //!< Level for debug information
 
/** @name Constructors and destructors */
/*        ---------------------------- */
//@{
   //! Default creation constructor (assign default values)
   Euler3D_ThermallyPerfect_cState(void): 
    rho(DENSITY_STDATM), E(PRESSURE_STDATM/(rho*(0.4))) {
     rhospec_null(); rhov.zero(); set_initial_values(); 
   }   

   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &value) : 
    rho(value), E(value) {
      rhov.x = value; rhov.y = value; rhov.z = value;
      rhospec_null(); set_initial_values(value); 
   }   
  
   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &V, 
                                   const double &En) :
     rho(d), rhov(V), E(En) {
      rhospec_null(); set_initial_values();
   } 

   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &V, 
                                   const double &En, 
                                   const double &frac) :
    rho(d), rhov(V), E(En) {
      rhospec_null(); set_initial_values(frac);
   } 

   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const Vector3D &dV, 
                                   const double &En, 
                                   const Species *rhomfrac) :
    rho(d), rhov(dV), E(En) {
      rhospec_null(); set_initial_values(rhomfrac);
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz,
                                   const double &En) :
    rho(d), E(En) {
      rhov.x=vx; rhov.y=vy; rhov.z=vz; 
      rhospec_null(); set_initial_values(); 
   }
  
   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &En,	
                                   const Species *rhomfrac) : 
    rho(d), E(En) { 
      rhov.x=vx; rhov.y=vy; rhov.z = vz;
      rhospec_null(); set_initial_values(rhomfrac); 
   }

   //! Assignment constructor
   Euler3D_ThermallyPerfect_cState(const double &d, 
                                   const double &vx, 
                                   const double &vy, 
                                   const double &vz, 
                                   const double &En,	
                                   const double &rhomfrac) : 
    rho(d), E(En) { 
      rhov.x=vx; rhov.y=vy; rhov.z = vz;
      rhospec_null(); set_initial_values(rhomfrac); 
   }
  
   //! Copy constructor (this is needed for the operator overload returns)
   Euler3D_ThermallyPerfect_cState(const Euler3D_ThermallyPerfect_cState &U) { 
      rhospec_null(); rho = DENSITY_STDATM; set_initial_values(); Copy(U); 
   }  

   //! Default destructor
   ~Euler3D_ThermallyPerfect_cState(void){ 
     Deallocate();
   }
//@}

/** @name Some useful operators */
/*        --------------------- */
//@{
   //! Allocate memory for species data
   void Allocate(void) {
#ifdef STATIC_NUMBER_OF_SPECIES
      if (STATIC_NUMBER_OF_SPECIES < ns) {
        cerr << "\n WARNING USING STATIC ARRAYS FOR SPECIES DATA WITH STATIC_NUMBER_OF_SPECIES = " 
             << STATIC_NUMBER_OF_SPECIES << "." 
             << "\n HOWEVER THERE IS A REQUEST FOR FOR ns = " << ns << " SPECIES.\n"; 
        exit(1); 
      }
#else
      if (rhospec != NULL) {
         delete[] rhospec;
         rhospec = NULL;
      } /* endif */
      rhospec = new Species[ns];
#endif
   }

   //! Allocate memory for species data
   void Allocate(const int &n) {
#ifdef STATIC_NUMBER_OF_SPECIES
      if (STATIC_NUMBER_OF_SPECIES < n) {
        cerr << "\n WARNING USING STATIC ARRAYS FOR SPECIES DATA WITH STATIC_NUMBER_OF_SPECIES = " 
             << STATIC_NUMBER_OF_SPECIES << "." 
             << "\n HOWEVER THERE IS A REQUEST FOR FOR n = " << n << " SPECIES.\n"; 
        exit(1); 
      }
#else
      if (rhospec != NULL) {
         delete[] rhospec;
         rhospec = NULL;
      } /* endif */
      ns = n;
      rhospec = new Species[ns];
#endif
   }

   //! Deallocate memory for species data
   void Deallocate(void) {
#ifndef STATIC_NUMBER_OF_SPECIES 
      if (rhospec != NULL) delete[] rhospec;
      rhospec = NULL;
#endif
   }

   //! Allocate static memory for species data
   void Allocate_static(void) {
      if (_temp_values != NULL) {
         delete[] _temp_values;
         _temp_values = NULL;
      } /* endif */
      _temp_values = new double[ns];
      if (_diff_coeff != NULL) {
         delete[] _diff_coeff;
         _diff_coeff = NULL;
      } /* endif */
      _diff_coeff = new double[ns];
      if (specdata != NULL) {
         delete[] specdata;
         specdata = NULL;
      } /* endif */
      specdata = new NASARP1311data[ns]; 
      if (Schmidt != NULL) {
         delete[] Schmidt; 
         Schmidt = NULL;
      } /* endif */
      Schmidt = new double[ns];
   }

   //! Allocate static memory for species data
   void Allocate_static(const int &n) {
      ns = n;
      if (_temp_values != NULL) {
         delete[] _temp_values;
         _temp_values = NULL;
      } /* endif */
      _temp_values = new double[ns];
      if (_diff_coeff != NULL) {
         delete[] _diff_coeff;
         _diff_coeff = NULL;
      } /* endif */
      _diff_coeff = new double[ns];
      if (specdata != NULL) {
         delete[] specdata;
         specdata = NULL;
      } /* endif */
      specdata = new NASARP1311data[ns]; 
      if (Schmidt != NULL) {
         delete[] Schmidt; 
         Schmidt = NULL;
      } /* endif */
      Schmidt = new double[ns];
   }

   //! Deallocate static memory for species data
   void Deallocate_static(void) {
      if (_temp_values != NULL) delete[] _temp_values;
      _temp_values = NULL;
      if (_diff_coeff != NULL) delete[] _diff_coeff;
      _diff_coeff = NULL;
      if (specdata != NULL) delete[] specdata;
      specdata = NULL;
      if (Schmidt != NULL)  delete[] Schmidt;
      Schmidt = NULL;
   }

   //! Set species point to null
   void rhospec_null(void) {
#ifndef STATIC_NUMBER_OF_SPECIES 
      rhospec = NULL;
#endif
   }

   //! Assign the species data (needs to be called only once as it is static)
   void set_species_data(const int &n,
                         const string *S,
                         const char *PATH,
                         const int &debug, 
                         const double &Mr, 
                         const double* Sc,
                         const int &trans_data);

   //! Set initial values for species data
   void set_initial_values(void);
   //! Set initial values for species data
   void set_initial_values(const double &value);
   //! Set initial values for species data
   void set_initial_values(const double *rhomfrac);
   //! Set initial values for species data
   void set_initial_values(const Species *rhomfrac);
     
   //! Copy solution state (cheaper than = operator)
   void Copy(const Euler3D_ThermallyPerfect_cState &U);

   //! Assigns a vacuum solution state
   void Vacuum(void){
      rho=ZERO; rhov.x=ZERO; rhov.y=ZERO; rhov.z=ZERO; E=ZERO;
      for (int i = 0; i < ns; i++) {
         rhospec[i].Vacuum();
      } /* endfor */
   }  

   //! Set the lower bound for valid temperature ranges for thermodynamic data
   void Temp_low_range(void);     
   //! Set the upper bound for valid temperature ranges for thermodynamic data 
   void Temp_high_range(void); 

   //! Check for negative species mass fractions
   bool negative_speccheck(void);

   //! Check for physical validity of the solution vector
   bool Realizable_Solution_Check(void);
//@}

/** @name Thermodynamic and other state functions */
/*        --------------------------------------- */
//@{
   //! Mixture gas constant
   double Rtot(void) const; 

   //! Mixture heat of formation
   double HeatofFormation(void) const;

   //! Mixture specific heat at constant pressure
   double Cp(void) const;  

   //! Mixture specific heat at constant pressure
   double Cv(void) const;  

   //! Mixture specific heat ratio
   double g(void) const;   

   //! Estimate for mixture specific heat ratio
   double gamma_guess(void) const;  

   //! Mixture absolute (sensible+chemical) internal energy
   double e(void) const;            

   //! Mixture sensible internal energy
   double es(void) const;           

   //! Mixture absolute specific enthalpy  
   double h(void) const; 
   //! Mixture absolute specific enthalpy  
   double h(const double &T) const;  

   //! Mixture sensible specific enthalpy
   double hs(void) const;
   //! Mixture sensible specific enthalpy
   double hs(const double &T) const; 

   //! Derivative of specific enthalpy wrt temperature, dh/dT
   double hprime(const double &T) const; 

   //! Mixture velocity
   Vector3D v(void) const;   

   //! Mixture pressure
   double p(void) const;    

   //! Mixture temperature
   double T(void) const;      

   //! Mixture sound speed
   double a(void) const;

   //! Mixture Mach number
   double M(void) const;   

   //! Determine the sum of species mass fractions
   double sum_species(void) const;
//@}
 
/** @name Primitive solution state */ 
/*        ------------------------ */
//@{
   //! Returns primitive solution state
   Euler3D_ThermallyPerfect_pState W(void);
   //! Returns primitive solution state
   Euler3D_ThermallyPerfect_pState W(void) const;
   //! Returns primitive solution state
   Euler3D_ThermallyPerfect_pState W(const Euler3D_ThermallyPerfect_cState &U) const;
//@}
 
/** @name Numerical Flux Functions */
/*        ------------------------ */
//@{
   //! Returns rotated conserved state aligned with x-axis in the norm_dir
   Euler3D_ThermallyPerfect_cState Rotate(const Vector3D &norm_dir) const;

   //! Returns un-rotated conserved state aligned with x-axis of the global problem
   Euler3D_ThermallyPerfect_cState RotateBack(const Vector3D &norm_dir) const;
//@}

/** @name Operators */
/*        --------- */
//@{
   //! Index operator
   double &operator[](int index);
   //! Index operator
   const double &operator[](int index) const;
 
   //! Binary addition operator
   Euler3D_ThermallyPerfect_cState operator +(const Euler3D_ThermallyPerfect_cState &U) const;

   //! Binary subtraction operator
   Euler3D_ThermallyPerfect_cState operator -(const Euler3D_ThermallyPerfect_cState &U) const;

   //! Binary multiplication operator
   Euler3D_ThermallyPerfect_cState operator *(const double &a) const;

   //! Binary multiplication operator
   friend Euler3D_ThermallyPerfect_cState operator *(const double &a, 
                                                     const Euler3D_ThermallyPerfect_cState &U);

   //! Binary multiplication operator
   double operator *(const Euler3D_ThermallyPerfect_cState &U) const;

   //! Binary division operator
   Euler3D_ThermallyPerfect_cState operator /(const double &a) const;

   //! Binary vector product operator
   Euler3D_ThermallyPerfect_cState operator ^(const Euler3D_ThermallyPerfect_cState &U) const;

   //! Assignment operator
   Euler3D_ThermallyPerfect_cState& operator =(const Euler3D_ThermallyPerfect_cState &U); 

   //! Shortcut addition operator
   Euler3D_ThermallyPerfect_cState& operator +=(const Euler3D_ThermallyPerfect_cState &U);

   //! Shortcut subtraction operator
   Euler3D_ThermallyPerfect_cState& operator -=(const Euler3D_ThermallyPerfect_cState &U);
      
   //! Unary subtraction operator
   friend Euler3D_ThermallyPerfect_cState operator -(const Euler3D_ThermallyPerfect_cState &U);
  
   //! Equal relational operator
   friend int operator ==(const Euler3D_ThermallyPerfect_cState &U1, 
                          const Euler3D_ThermallyPerfect_cState &U2);

   //! Not equal relational operator
   friend int operator !=(const Euler3D_ThermallyPerfect_cState &U1, 
                          const Euler3D_ThermallyPerfect_cState &U2);
  
   //! Output stream operator
   friend ostream& operator << (ostream &out_file, 
                                const Euler3D_ThermallyPerfect_cState &U);

   //! Input stream operator
   friend istream& operator >> (istream &in_file,  
                                Euler3D_ThermallyPerfect_cState &U);
//@}   
};

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState member functions                                    *
 ***************************************************************************************/

/***************************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Index operators.                                 *
 ***************************************************************************************/
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

/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Binary arithmetic operators.        *
 **************************************************************************/
//------------------ Addition ------------------------//
inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
operator +(const Euler3D_ThermallyPerfect_pState &W) const { 
   Euler3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
   Temp.Copy(*this);
   Temp += W;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
operator -(const Euler3D_ThermallyPerfect_pState &W) const {
   Euler3D_ThermallyPerfect_pState Temp(W.rho,W.v,W.p);
   Temp.Copy(*this);
   Temp -= W;
   return Temp;
}

//---------------- Scalar Multiplication ------------------//
inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
operator *(const double &a) const {
   Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.v = v*a; Temp.p = p*a;
   for (int i=0; i<ns; i++) {
    Temp.spec[i] = spec[i]*a;
   } /* endif */ 
   return(Temp);
}

inline Euler3D_ThermallyPerfect_pState operator *(const double &a, 
                                                  const Euler3D_ThermallyPerfect_pState &W) {
   Euler3D_ThermallyPerfect_pState Temp;
   //Temp.Copy(W);
   Temp.rho = W.rho*a;  Temp.v = W.v*a; Temp.p = W.p*a;
   for (int i=0; i<W.ns; i++) {
      Temp.spec[i] = W.spec[i]*a;
   } /* endfor */
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
operator /(const double &a) const {
   Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.v = v/a; Temp.p = p/a; 
   for (int i = 0; i < ns; i++) {
     Temp.spec[i] = spec[i]/a; 
   } /* endfor */
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double Euler3D_ThermallyPerfect_pState::
operator *(const Euler3D_ThermallyPerfect_pState &W) const {
   double sum=0.0;
   for (int i = 0; i < ns; i++) {
      sum += spec[i]*W.spec[i];
   } /* endfor */
   return (rho*W.rho + v*W.v + p*W.p + sum);
}

//----------- Solution state product operator ------------//
inline Euler3D_ThermallyPerfect_pState Euler3D_ThermallyPerfect_pState::
operator ^(const Euler3D_ThermallyPerfect_pState &W) const {
   Euler3D_ThermallyPerfect_pState Temp(rho,v,p);
   Temp.Copy(*this);
   Temp.rho = rho*W.rho;
   Temp.v.x = v.x*W.v.x;
   Temp.v.y = v.y*W.v.y;
   Temp.v.z = v.z*W.v.z;
   Temp.p = p*W.p;
   for (int i=0; i<ns; i++) {
      Temp.spec[i] = spec[i]*W.spec[i];
   } /* endif */ 
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::
operator =(const Euler3D_ThermallyPerfect_pState &W) {
   //self assignment protection
   if (this != &W) {   
      //copy assignment
      rho = W.rho;
      v = W.v; 
      p = W.p; 
      for (int i=0; i < ns; i++) {
	spec[i] = W.spec[i];
      } /* endfor */
   } /* endif */
   return (*this);
}

/*************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Shortcut arithmetic operators.     *
 *************************************************************************/
inline Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::
operator +=(const Euler3D_ThermallyPerfect_pState &W) {
   rho += W.rho;
   v += W.v; 
   p += W.p; 
   for (int i = 0; i < ns; i++) {
      spec[i].c += W.spec[i].c;
   } /* endfor */
   return (*this);
}

inline Euler3D_ThermallyPerfect_pState& Euler3D_ThermallyPerfect_pState::
operator -=(const Euler3D_ThermallyPerfect_pState &W) {
   rho -= W.rho;
   v -= W.v;
   p -= W.p;
   for (int i = 0; i < ns; i++) {
      spec[i].c -= W.spec[i].c;
   } /* endfor */ 
   return (*this); 
}

/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Unary arithmetic operators.         *
 **************************************************************************/
inline Euler3D_ThermallyPerfect_pState operator -(const Euler3D_ThermallyPerfect_pState &W) {
   Species *spt= new Species[W.ns];
   for (int i = 0; i < W.ns; i++) {
      spt[i] = -W.spec[i]; 
   } /* endfor */ 
   Euler3D_ThermallyPerfect_pState Temp(-W.rho, -W.v, -W.p, spt);
   delete[] spt;
   return(Temp);
}

/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Relational operators.               *
 **************************************************************************/
inline int operator ==(const Euler3D_ThermallyPerfect_pState &W1, 
                       const Euler3D_ThermallyPerfect_pState &W2) {
   bool Temp;
   for (int i = 0; i < W1.ns; i++) {
      if (W1.spec[i] == W2.spec[i] ){
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */  
   } /* endfor */
   return (W1.rho == W2.rho && W1.v == W2.v && W1.p == W2.p && Temp == true);
}

inline int operator !=(const Euler3D_ThermallyPerfect_pState &W1, 
                       const Euler3D_ThermallyPerfect_pState &W2) {
    bool Temp = true;
    for (int i=0; i<W1.ns; i++){
       if (W1.spec[i] != W2.spec[i]) {
         Temp = false;
         break;
       } /* endif */
    } /* endfor */
    return (W1.rho != W2.rho || W1.v != W2.v || W1.p != W2.p || Temp != true);
}

/**************************************************************************
 * Euler3D_ThermallyPerfect_pState -- Input-output operators.             *
 **************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Euler3D_ThermallyPerfect_pState &W) {
   out_file.precision(10);
   out_file.setf(ios::scientific);
   out_file << " " << W.rho  << " " << W.v.x << " " << W.v.y 
            << " " << W.v.z << " " << W.p;
   for (int i=0; i<W.ns; i++) {
      out_file<<" "<<W.spec[i];
   } /* endfor */
   out_file.unsetf(ios::scientific);
   return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             Euler3D_ThermallyPerfect_pState &W) {
   in_file.setf(ios::skipws);
   in_file >> W.rho >> W.v.x >> W.v.y >> W.v.z >>W.p;
   //W.set_initial_values();
   for (int i = 0; i < W.ns; i++) {
      in_file>>W.spec[i];
   } /* endfor */
   in_file.unsetf(ios::skipws);
   return (in_file);
}

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState member functions                                    *
 ***************************************************************************************/

/***************************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Index operators.                                 *
 ***************************************************************************************/
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

inline const double& Euler3D_ThermallyPerfect_cState::operator[](int index) const {
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

/**************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Binary arithmetic operators.        *
 **************************************************************************/
//----------------- Addition -----------------------------//
inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator +(
   const Euler3D_ThermallyPerfect_cState &U) const{ 
   Euler3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
   Temp.Copy(*this);
   Temp += U;
   return Temp;
}

//------------------ Subtraction ------------------------//
inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::operator -(
   const Euler3D_ThermallyPerfect_cState &U) const{
   Euler3D_ThermallyPerfect_cState Temp(U.rho,U.rhov,U.E);
   Temp.Copy(*this);
   Temp -= U;
   return Temp;
}


//---------------- Scalar Multiplication ------------------//
inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::
operator *(const double &a) const {
   Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho*a;  Temp.rhov = rhov*a; Temp.E = E*a;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rhospec[i]*a;
   } /* endfor */
   return(Temp);
}

inline Euler3D_ThermallyPerfect_cState operator *(const double &a, 
                                                  const Euler3D_ThermallyPerfect_cState &U) {
   Euler3D_ThermallyPerfect_cState Temp;
   Temp.rho = U.rho*a;  Temp.rhov = U.rhov*a; Temp.E = U.E*a;
   for (int i=0; i<U.ns; i++){
      Temp.rhospec[i] = U.rhospec[i]*a;
   } /* endfor */
   return(Temp);
}

//--------------- Scalar Division ------------------------//
inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::
operator /(const double &a) const {
   Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho/a; Temp.rhov = rhov/a; Temp.E = E/a;
   for (int i = 0; i < ns; i++) {
      Temp.rhospec[i] = rhospec[i]/a; 
   } /* endfor */
   return(Temp);
}

//----------------- Inner Product ------------------------//
inline double Euler3D_ThermallyPerfect_cState::
operator *(const Euler3D_ThermallyPerfect_cState &U) const {
   double sum=0.0;
   for (int i = 0; i < ns; i++) {
      sum += rhospec[i]*U.rhospec[i];
   } /* endfor */
   return (rho*U.rho + rhov*U.rhov + E*U.E + sum);
}

//----------- Solution state product operator ------------//
inline Euler3D_ThermallyPerfect_cState Euler3D_ThermallyPerfect_cState::
operator ^(const Euler3D_ThermallyPerfect_cState &U) const {
   Euler3D_ThermallyPerfect_cState Temp(rho,rhov,E);
   Temp.Copy(*this);
   Temp.rho = rho*U.rho;
   Temp.rhov.x = rhov.x*U.rhov.x;
   Temp.rhov.y = rhov.y*U.rhov.y;
   Temp.rhov.z = rhov.z*U.rhov.z;
   Temp.E = E*U.E;
   for(int i=0; i<ns; i++){
      Temp.rhospec[i] = rhospec[i]*U.rhospec[i];
   }  
   return(Temp);
}

//----------------- Assignment ----------------------------//
inline Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::
operator =(const Euler3D_ThermallyPerfect_cState &U) {
   //self assignment protection
   if (this != &U) {   
      //copy assignment
      rho = U.rho;
      rhov = U.rhov; 
      E = U.E; 
      for (int i = 0; i < ns; i++) {
         rhospec[i] = U.rhospec[i];
      } /* endfor */
   } /* endif */
   return (*this);
}

/*************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Shortcut arithmetic operators.     *
 *************************************************************************/
inline Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::
operator +=(const Euler3D_ThermallyPerfect_cState &U) {
   rho += U.rho;
   rhov += U.rhov; 
   E += U.E;
   for(int i = 0; i < ns; i++) {
      rhospec[i].c += U.rhospec[i].c;
   } /* endfor */
   return (*this);
}

inline Euler3D_ThermallyPerfect_cState& Euler3D_ThermallyPerfect_cState::
operator -=(const Euler3D_ThermallyPerfect_cState &U) {
   rho -= U.rho;
   rhov -= U.rhov;
   E -= U.E;
   for (int i = 0; i < ns; i++) {
     rhospec[i].c -= U.rhospec[i].c;
   } /* endfor */
   return (*this); 
}

/*************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Unary arithmetic operators.        *
 *************************************************************************/
inline Euler3D_ThermallyPerfect_cState operator -(const Euler3D_ThermallyPerfect_cState &U) {
   Species *spt= new Species[U.ns];
   for(int i=0; i<U.ns; i++){
      spt[i] = -U.rhospec[i]; 
   }  
   Euler3D_ThermallyPerfect_cState Temp(-U.rho,-U.rhov,-U.E, spt);
   delete[] spt;
   return(Temp);
}

/*************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Relational operators.              *
 *************************************************************************/
inline int operator ==(const Euler3D_ThermallyPerfect_cState &U1, 
                       const Euler3D_ThermallyPerfect_cState &U2) {
   bool Temp;
   for (int i = 0; i < U1.ns; i++) {
      if( U1.rhospec[i] == U2.rhospec[i] ) {
         Temp = true;
      } else {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho == U2.rho && U1.rhov == U2.rhov && U1.E == U2.E && Temp == true);
}

inline int operator !=(const Euler3D_ThermallyPerfect_cState &U1, 
                       const Euler3D_ThermallyPerfect_cState &U2) {
   bool Temp = true;
   for (int i = 0; i < U1.ns; i++) {
      if (U1.rhospec[i] != U2.rhospec[i]) {
         Temp = false;
         break;
      } /* endif */
   } /* endfor */
   return (U1.rho != U2.rho || U1.rhov != U2.rhov || U1.E != U2.E || Temp != true);
}

/*************************************************************************
 * Euler3D_ThermallyPerfect_cState -- Input-output operators.            *
 *************************************************************************/
inline ostream &operator << (ostream &out_file, 
                             const Euler3D_ThermallyPerfect_cState &U) {
  //out_file.precision(20);
  out_file.setf(ios::scientific);
  out_file << " " << U.rho  << " " << U.rhov.x << " " << U.rhov.y << " " 
           << U.rhov.z << " "<< U.E;
  for(int i = 0; i < U.ns; i++) {
     out_file<<" "<<U.rhospec[i];
  } /* endfor */
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, 
                             Euler3D_ThermallyPerfect_cState &U) {
   in_file.setf(ios::skipws);
   in_file >> U.rho >> U.rhov.x >> U.rhov.y >> U.rhov.z >> U.E;
  //U.set_initial_values();
  for( int i = 0; i < U.ns; i++) {
    in_file>>U.rhospec[i]; 
  } /* endfor */
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif // _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
