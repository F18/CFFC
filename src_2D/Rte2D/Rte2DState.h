/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Rte2DState.h                                               **
 **                                                                  **
 ** Description: The radiation state class contains the properties   **
 **              of a gray, absorbing, emitting, anisotropically     **
 **              scattering medium.  Additionally, it contains the   **
 **              state (intensity) of a beam of light passing        **
 **              through the medium.  This file defines the state    ** 
 **              class.                                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _RTE2D_STATE_INCLUDED
#define _RTE2D_STATE_INCLUDED 

/********************************************************
 * Class Declaration                                    *
 ********************************************************/
class Rte2D_State;

/********************************************************
 * Include required C++ libraries                       *
 ********************************************************/
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/********************************************************
 * Required CFDkit+caboodle header files                *
 ********************************************************/

#include "../Math/Math.h"
#include "../Math/Matrix.h"
#include "../Math/Vector2D.h"
#include "../Math/Vector3D.h"
#include "../Math/Simpson.h"
#include "../CFD/CFD.h"
#include "../Physics/GasConstants.h"
#include "SNBCK.h"
#include "Planck.h"

/********************************************************
 * Necessary Rte2D Specific Constants                   *
 ********************************************************/

// Absorbsion model type 
enum Gas_Models { RTE2D_ABSORB_GRAY,
                  RTE2D_ABSORB_SNBCK };

// RTE discretization Type (DOM or FVM)
enum RTE_Discretization { RTE2D_SOLVER_FVM,
			  RTE2D_SOLVER_DOM };

// Scatter Phase Functions
enum Scatter_Models { RTE2D_SCATTER_ISO,  // isotropic scattering
		      RTE2D_SCATTER_F1,   // Kim and Lee (1988)   
		      RTE2D_SCATTER_F2,   // Kim and Lee (1988)
		      RTE2D_SCATTER_F3,   // Kim and Lee (1988) 
		      RTE2D_SCATTER_B1,   // Kim and Lee (1988)
		      RTE2D_SCATTER_B2 }; // Kim and Lee (1988)

// DOM quadrature Type 
enum DOM_Qauad { RTE2D_DOM_S2,   // S2 NONSYMMETRIC from LATHROP and CARLSON
		 RTE2D_DOM_S4,   // S4 SYMMETRIC from LATHROP and CARLSON
		 RTE2D_DOM_S6,   // S6 SYMMETRIC from LATHROP and CARLSON
		 RTE2D_DOM_S8,   // S8 SYMMETRIC from LATHROP and CARLSON
		 RTE2D_DOM_S12,  // S12 SYMMETRIC from LATHROP and CARLSON
		 RTE2D_DOM_T3 }; // T3 SYMMETRIC from Truelove


/********************************************************
 * STRUCTS REQUIRED FOR INTEGRATION                     *
 ********************************************************/

// struct needed to integrate the phase function over the solid angle
// contains information for legendre polynomials
struct legendre_param {
  double *An; // the expansion coefficient array
  int Mn;     // degree of Legendre polynomial 
};

// struct needed to integrate the exact solution for radiative
// heat transfer in a cylindrical enclosure
struct exact_cyl_param {
  double z;       // the non-dimensional axial distance
  double r;       // the non-dimensional radial distance 
  double c;       // the half length divided by the outer radius
  double kappa;   // non-dimensional absorbsion coefficient
  int term_flag;  // a flag for which parameter we are computing
  int coord_flag; // a flag for which coordinate this is for
};

// struct needed to integrate the exact solution for radiative
// heat transfer in a rectangular enclosure
struct exact_rect_param {
  double x;       // the dimensional x location
  double y;       // the dimensional y location 
  double a1;      // west wall location
  double a2;      // east wall location
  double b1;      // south wall location
  double b2;      // north wall location
  double kappa;   // absorbsion coefficient
  int term_flag;  // a flag for which parameter we are computing
  int coord_flag; // a flag for which coordinate this is for
};


/*!
 * Class: Rte2D_State
 *
 * @brief Solution state class definition for am infinitesmally thin 
 *        pencil of rays through a non-gray mixture of absorbing,
 *        emitting, anisotropically scattering gases.
 *
 * Solution state class definition for am infinitesmally thin 
 * pencil of rays through a non-gray mixture of absorbing, emitting, 
 * anisotropically scattering gases.  
 *
 * \verbatim
 * Member functions
 *     I        -- Return array of spectral, directional intensities
 *     kappa    -- Return array of spectral absorbsion coefficient
 *     sigma    -- Return array of spectral scattering coefficient
 *     Ib       -- Return array of spectral blackbody intensity.
 *     Npolar   -- Return the number of polar directions.
 *     Nazim    -- Return array of the number of azimuthal directions in each polar direction
 *     Nbands   -- Return the number of spectral frequency bands.
 *     Ntot     -- Total number of directions times total number of bands.
 *     Index    -- Index array relating direction, and frequency indices to unfolded array of intensities.
 *     NUM_VAR_RTE2D -- Returns the number of variables in the solution state.
 *     mu       -- Returns the (x,r)-direction cosine.
 *     eta      -- Returns the (y,azimuthal)-direction cosine.
 *     xi       -- Returns the (z)-direction cosine.
 *     omega    -- Return the array of discrete control angle element sizes
 *     theta    -- Return the arraoy of polar angle grid points
 *     psi      -- Return the array of azimuthal angle grid points
 *     delta_theta     -- Return the array of polar angle grid points
 *     delta_psi       -- Return the array of azimuthal angle grid points
 *     Phi      -- Return the array of scattering phase function
 *     Symmetry_Factor -- Return the solid angle range symmetry factor
 *     RTE_Type -- Return the flag for DOM or FVM solver type
 *     I_half   -- Returns storage array for intensity in special directions
 *     alpha    -- Returns angular redistribution coefficient (for DOM)
 *     U        -- Return conserved solution state.
 *     Fn       -- Return n-direction solution flux.
 *     dFndU    -- Return n-direction flux Jacobian.
 *     S        -- Return regular source term.
 *     dSdU     -- Return regular source term jacobian.
 *     Sa       -- Return axisymmetric source terms (Time-marching case).
 *     Sa_FVM   -- Return axisymmetric source terms (FVM Space-marching case).
 *     Sa_DOM   -- Return axisymmetric source terms (DOM Space-marching case).
 *     dSadU    -- Return axisymmetric source term jacobian (Time-marching case).
 *     dSadU_FVM-- Return axisymmetric source term jacobian (FVM Space-marching case).
 *     dSadU_DOM-- Return axisymmetric source term jacobian (DOM Space-marching case).
 *     In       -- Return intensity at spectral/directional location.
 *     beta     -- Return extinction coefficient.
 *     G        -- Return total directional integrated radiation.
 *     q        -- Return heat flux vector.
 *
 * Member operators
 *      U -- a conserved solution state
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

class Rte2D_State {
 private:
 public:

  //@{ @name Conserved variables and associated constants:
  double* kappa;              //!< absorbsion coefficient [m^-1]
  double* sigma;              //!< scattering coefficient [m^-1]
  double* Ib;                 //!< blackbody intentsity [W/m^2]
  double* I;                  //!< directional intentsity [W/m^2]
  //@}

  //@{ @name DOM space marching parameters:
  double** I_half;            //!< storage for intensity in special directions.  See Carlson and Lathrop (1968)
  double** alpha;             //!< angular redistribution coefficient (for DOM)
  //@}

  //@{ @name Static index variables:
  static int  Npolar;         //!< number of polar directions
  static int* Nazim;          //!< number of azimuthal directions
  static int  Nband;          //!< the total number of frequency bands
  static int  Ntot;           //!< total number of intensities
  static int*** Index;        //!< indexing array
  static int NUM_VAR_RTE2D;   //!< total number of Rte2D variables
  //@}

  //@{ @name Static variables related to angular discretization:
  static double** mu;         //!< (x,r)-direction cosine
  static double** eta;        //!< (y,azimuthal)-direction cosine
  static double** xi;         //!< (z)-direction cosine
  static double** omega;      //!< discrete control angle element sizes
  static double* theta;       //!< polar angle grid points
  static double** psi;        //!< azimuthal angle grid points
  static double* delta_theta; //!< polar angle grid points
  static double** delta_psi;  //!< azimuthal angle grid points

  //@{ @name Miscillaneous static variables:
  static double***** Phi;     //!< scattering phase function
  static int Symmetry_Factor; //!< solid angle range symmetry factor
  static SNBCK* SNBCKdata;    //!< statistical narrow band model 
  //@}


  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Rte2D_State() : kappa(NULL), sigma(NULL), I(NULL), Ib(NULL), 
                  alpha(NULL), I_half(NULL)
    { Allocate(); }

  //! Copy constructor.
  Rte2D_State( const Rte2D_State &U ) : kappa(NULL), 
                                        sigma(NULL),			     
                                        I(NULL), 
                                        Ib(NULL),
                                        alpha(NULL), 
                                        I_half(NULL)
    { Allocate(); if( this != &U) Copy(U); }
  
  //! Destructor.
  ~Rte2D_State() { Deallocate(); }
  //@}


  //@{ @name Useful operators.
  //! Copy solution state operator.
  void Copy( const Rte2D_State &U ) {
    for ( int i=0; i<Nband; i++ ) { 
      Ib[i] = U.Ib[i];
      kappa[i] = U.kappa[i];
      sigma[i] = U.sigma[i];
    }
    for ( int i=0; i<Ntot; i++ ) I[i] = U.I[i];
  }

  //! Copy non solution state operator.
  void Copy_NonSol( const Rte2D_State &U ) {
    for ( int i=0; i<Nband; i++ )   { 
      Ib[i] = U.Ib[i];
      kappa[i] = U.kappa[i];
      sigma[i] = U.sigma[i];
    }
  }

  //! Vacuum operator.
  void Vacuum() {
    for ( int i=0; i<Ntot; i++ ) I[i] = ZERO;
    zero_non_sol();
  }

  //! Zero solution operator.
  void Zero() {
    for ( int i=0; i<Ntot; i++ ) I[i] = ZERO;
    zero_non_sol();
  }

  //! Zero non-solution operator.
  void zero_non_sol() {
    for(int i=0; i<Nband; i++) {
      kappa[i] = 0;
      sigma[i] = 0;
      Ib[i] = 0;  
    }
  }  

  //! Check for unphysical state properties.
  int Unphysical_Properties(void) const {
    for ( int n=0 ; n<Ntot ; n++ ) if ( I[n]<ZERO || I[n]!=I[n] ) return 1;
    return 0;
  }
  //@}


  //@{ @name Initialization functions
  //! Set the spectral/directional intensity to a constant value.
  void SetIntensity(double val);

  //! Set the spectral absorbsion coefficient.
  void SetAbsorption(double val);

  //! Set the spectral scattering coefficient.
  void SetScattering(double val);

  //! Set the spectral blackbody intensity.
  void SetBlackbody(double val);

  //! Compute DOM angular redistribution coefficients
  void SetupART_DOM( const Vector2D nfaceE, const double AfaceE,
		     const Vector2D nfaceW, const double AfaceW,
		     const Vector2D nfaceN, const double AfaceN,
		     const Vector2D nfaceS, const double AfaceS );
  //@}

  //@{ @name Functions to setup the absoorption coefficient model.
  //! Initialize model type.
  static void SetupAbsorb( const SNBCK_Input_Parameters &IP, 
			   const char* CFFC_PATH );
  //@}


  //@{ @name Functions to set the number of directions and compute the cosines.
  //! For FVM
  static void SetDirsFVM(const int NumPolarDirs, const int NumAzimDirs, 
			 const int Axisymmetric);
  //! For DOM
  static void SetDirsDOM(const int Quad_Type, const int Axisymmetric,
			 const char *CFFC_PATH );
  //@}


  //@{ @name Functions to setup the phase function.
  //! For FVM
  static void SetupPhaseFVM( const int type );
  //! For DOM
  static void SetupPhaseDOM( const int type );
  //@}


  //@{ @name Allocators and deallocators
  //! memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();

  //! memory allocation / deallocation for the directional arrays
  static void AllocateDirs( );
  static void DeallocateDirs();

  //! memory allocation / deallocation for the cosine arrays
  static void AllocateCosines( const int RTE_Type );
  static void DeallocateCosines();

  //! memory allocation / deallocation for the SNBCKdata object
  static void AllocateSNBCK();
  static void DeallocateSNBCK();

  //! deallocate all static objects
  static void DeallocateStatic() 
    { DeallocateCosines(); DeallocateDirs(); DeallocateSNBCK(); }

  //@}


  //@{ @name Multigrid operators.
  //! Copy variables solved by multigrid only.
  void Copy_Multigrid_State_Variables(const Rte2D_State &Ufine) {
    Copy(Ufine);
  }

  //! Zero variables not-solved by multigrid.
  void Zero_Non_Multigrid_State_Variables(void) { }
  //@}


  //@{ @name State functions.
  //! access intensity using 2D indexes 
  double& In( const int v, const int m, const int l ) const;
  double& In( const int v, const int m, const int l );

  //! Return extinction coefficient
  double beta(const int v) const;

  //! Return directional integrated intensity [W/m^2]
  double G();        

  //! Return heat flux vector [W/m^2]
  Vector2D q();
  //@}

  //@{ @name Solution flux and Jacobian (n-direction).
  //! Flux - For time-marching.
  Rte2D_State Fn( const Vector2D &norm_dir ) const;
  //! Flux - For space-marching
  double Fn( const Vector2D &norm_dir, const int v, 
	     const int m, const int l ) const;
  //! Flux Jac. - For time-marching.
  friend void dFndU(DenseMatrix &dFdU, const Rte2D_State &U, 
		    const Vector2D &norm_dir);
  //! Flux Jac. - For space-marching
  double dFndU( const Vector2D &norm_dir, const int m, const int l ) const;
  //@}

  //@{ @name Source vector (regular terms).
  //! Source - For time-marching.
  Rte2D_State S(void) const;
  //! Source - For space-marching.
  double S(const int v, const int m, const int l) const;
  //! Source Jac. - For time-marching.
  void dSdU(DenseMatrix &dSdU) const;
  //! Source Jac. - For space-marching.
  double dSdU(const int v, const int m, const int l ) const;
  //@}

  //@{ @name Source vector (axisymmetric terms).
  //! Source - For time-marching.
  Rte2D_State Sa( const Rte2D_State &dUdpsi,
		  const Rte2D_State &phi_psi,
		  const int Axisymmetric) const;
  //! Source - For space-marching (FVM).
  double Sa_FVM(const int v, const int m, const int l, 
		const int Axisymmetric) const;
  //! Source - For space-marching (DOM).
  double Sa_DOM(const int v, const int m, const int l) const;
  //! Source Jac. - For time-marching.
  void dSadU(DenseMatrix &dSadU, const double Sp,
	     const int Axisymmetric) const;
  //! Source Jac. - For space-marching (FVM).
  double dSadU_FVM(const int v, const int m, const int l,
		   const int Axisymmetric) const;
  //! Source Jac. - For space-marching (DOM).
  double dSadU_DOM(const int v, const int m, const int l) const;
  //@}


  //@{ @name Index operator.
  double &operator[](int index) { return I[index-1]; }
  const double &operator[](int index) const { return I[index-1]; }
  //@}

  //@{ @name Binary arithmetic operators.
  Rte2D_State operator +(const Rte2D_State &U) const;
  Rte2D_State operator -(const Rte2D_State &U) const;
  Rte2D_State operator *(const double &a) const;
  double operator   *(const Rte2D_State &U) const;
  friend Rte2D_State operator *(const double &a, const Rte2D_State &U);
  Rte2D_State operator /(const double &a) const;
  Rte2D_State operator /(const Rte2D_State &U) const;
  Rte2D_State operator ^(const Rte2D_State &U) const;
  //@}

  //@{ @name Assignment Operator. 
  Rte2D_State& operator =(const Rte2D_State &U); 
  //@}

  //@{ @name Shortcut arithmetic operators.
  Rte2D_State& operator +=(const Rte2D_State &U);
  Rte2D_State& operator -=(const Rte2D_State &U);
  Rte2D_State &operator *=(const double &a);
  Rte2D_State &operator /=(const double &a);
  //@}
      
  //@{ @name Unary arithmetic operators.
  friend Rte2D_State operator -(const Rte2D_State &U);
  //@}
  
  //@{ @name Relational operators.
  friend int operator ==(const Rte2D_State &U1, const Rte2D_State &U2);
  friend int operator !=(const Rte2D_State &U1, const Rte2D_State &U2);
  friend int operator >(const Rte2D_State &U1, const double d);
  friend int operator <(const Rte2D_State &U1, const double d);
  //@}

  //@{ @name Input-output operators.
  friend ostream& operator << (ostream &out_file, const Rte2D_State &U);
  friend istream& operator >> (istream &in_file,  Rte2D_State &U);
 //@}

};

 /**************************************************************************
  ********************* RTE2D_STATE ACCESSORS ******************************
 ***************************************************************************/


/********************************************************
 * Access intensity array using 3D indexing.            *
 ********************************************************/
inline double& Rte2D_State :: In( const int v, const int m, const int l ) const
{ return I[ Index[v][m][l] ]; }

inline double& Rte2D_State :: In( const int v, const int m, const int l )
{ return I[ Index[v][m][l] ]; }


/********************************************************
 * Return extinction coefficient.                       *
 ********************************************************/
inline double Rte2D_State :: beta(const int v) const
{ return (kappa[v] + sigma[v]); }


/*************************************************************
 *Compute directional integrated intensity*
 *                   in [W/m^2]                              *
 *************************************************************/
inline double Rte2D_State :: G( void ) 
{
  double sum = 0;
  for ( int v=0; v<Nband; v++ ) 
    for ( int m=0; m<Npolar; m++ ) 
      for ( int l=0; l<Nazim[m]; l++ ) 
	sum += omega[m][l] * In(v,m,l);
  return sum;
}

/*************************************************************
 * Rte2D_State::q -- Compute radiative heat flux vector      *
 *                   in [W/m^2]                              *
 *************************************************************/
inline Vector2D Rte2D_State :: q( void )
{
  Vector2D Temp(ZERO);
  for ( int v=0; v<Nband; v++ ) 
    for ( int m=0; m<Npolar; m++ ) 
      for ( int l=0; l<Nazim[m]; l++ ) {
	Temp.x +=  mu[m][l] * In(v,m,l);
	Temp.y += eta[m][l] * In(v,m,l);
    }
  return Temp;  
}

 /**************************************************************************
  ********************* RTE2D_PSTATE CONSTRUCTORS **************************
 ***************************************************************************/

/********************************************************
 * Set initial values.                                  *
 ********************************************************/

inline void Rte2D_State :: SetIntensity(double val) {
  for ( int i=0; i<Ntot; i++ ) I[i] = val;
}

inline void Rte2D_State :: SetAbsorption(double val) {
  for ( int i=0; i<Nband; i++ ) kappa[i] = val;
}

inline void Rte2D_State :: SetScattering(double val) {
  for ( int i=0; i<Nband; i++ ) sigma[i] = val;
}

inline void Rte2D_State :: SetBlackbody(double val) {
  for ( int i=0; i<Nband; i++ ) Ib[i] = val;
}



/********************************************************
 * Array allocator and deallocator for intensity array. *
 ********************************************************/
inline void Rte2D_State :: Allocate()
{
  // deallocate first
  Deallocate();

  // create the jagged array
  if (Ntot>0)  { I = new double[Ntot]; }
  if (Nband>0) {
    Ib = new double[Nband];
    kappa = new double[Nband];
    sigma = new double[Nband];
  }

}


inline void Rte2D_State :: Deallocate()
{
  if (     I != NULL ) { delete[] I;           I = NULL; }
  if ( kappa != NULL ) { delete[] kappa;   kappa = NULL; }
  if ( sigma != NULL ) { delete[] sigma;   sigma = NULL; }
  if (    Ib != NULL ) { delete[] Ib;         Ib = NULL; }
  if ( alpha != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] alpha[i];
    delete[] alpha; 
    alpha = NULL; 
  }
}


/********************************************************
 * Array allocator and deallocator for directional      *
 * arrays.                                              *
 ********************************************************/
inline void Rte2D_State :: AllocateDirs()
{
  // deallocate first
  DeallocateDirs();

  // create the array
  if (Npolar>0)  Nazim = new int[Npolar]; 

}


inline void Rte2D_State :: DeallocateDirs()
{
  if ( Nazim != NULL ) { delete[] Nazim; Nazim = NULL; }
}


/********************************************************
 * Array allocator and deallocator for cosine arrays.   *
 ********************************************************/
inline void Rte2D_State :: AllocateCosines(const int RTE_Type)
{
  // declare
  int cnt;

  // create the array
  if (Nband>0 && Npolar>0 && Nazim!=NULL) {  

    mu          = new double*[Npolar];
    eta         = new double*[Npolar];
    xi          = new double*[Npolar];
    omega       = new double*[Npolar];
    if (RTE_Type==RTE2D_SOLVER_FVM) {
      psi         = new double*[Npolar];
      delta_psi   = new double*[Npolar];
      theta       = new double [Npolar+1];
      delta_theta = new double [Npolar];
    } /* endif */

    for (int i=0; i<Npolar; i++ ) {
      mu[i]        = new double[ Nazim[i] ];
      eta[i]       = new double[ Nazim[i] ];
      xi[i]        = new double[ Nazim[i] ];
      omega[i]     = new double[ Nazim[i] ];
      if (RTE_Type==RTE2D_SOLVER_FVM) {
	psi[i]       = new double[ Nazim[i]+1 ];
	delta_psi[i] = new double[ Nazim[i] ];
      } /* endif */
    } /* endfor */
    
    // set the indexe array
    cnt = 0;
    Index = new int**[Nband];
    for (int v=0; v<Nband; v++ ) {

      Index[v] = new int*[Npolar];
      for (int m=0; m<Npolar; m++ ) {

	Index[v][m] = new int[ Nazim[m] ];
	for (int l=0; l<Nazim[m]; l++ ) {
	  Index[v][m][l] = cnt;
	  cnt++;
	} /* endfor */
      } /*endfor*/
    }/* endfor */


    // allocate the phase function
    Phi = new double****[Nband];
    for (int v=0; v<Nband; v++ ) {

      Phi[v] = new double***[Npolar];
      for (int m=0; m<Npolar; m++ ) {

	Phi[v][m] = new double**[ Nazim[m] ];
	for (int l=0; l<Nazim[m]; l++ ) {

	  Phi[v][m][l] = new double*[Npolar];
	  for (int p=0; p<Npolar; p++ ) {
	    
	    Phi[v][m][l][p] = new double[ Nazim[p] ];	    

	  } /* endfor */
	} /* endfor */
      } /*endfor*/
    }/* endfor */
 
  } /* endif */

}


inline void Rte2D_State :: DeallocateCosines()
{
  if ( mu != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] mu[i];
    delete[] mu; 
    mu = NULL; 
  } /* endif */
  if ( eta != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] eta[i];
    delete[] eta; 
    eta = NULL; 
  } /* endif */
  if ( xi != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] xi[i];
    delete[] xi; 
    xi = NULL; 
  } /* endif */
  if ( omega != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] omega[i];
    delete[] omega; 
    omega = NULL; 
  } /* endif */
  if ( psi != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] psi[i];
    delete[] psi; 
    psi = NULL; 
  } /* endif */
  if ( delta_psi != NULL ) { 
    for (int i=0; i<Npolar; i++ ) delete[] delta_psi[i];
    delete[] delta_psi; 
    delta_psi = NULL; 
  } /* endif */
  if ( theta != NULL ) { 
    delete[] theta; 
    theta = NULL; 
  } /* endif */
  if ( delta_theta != NULL ) { 
    delete[] delta_theta; 
    delta_theta = NULL; 
  } /* endif */
  if ( Index != NULL ) { 
    for (int v=0; v<Nband; v++) {
      for (int m=0; m<Npolar; m++ ) delete[] Index[v][m];
      delete[] Index[v];
    } /* endfor */
    delete[] Index; 
    Index = NULL; 
  } /* endif */
  if ( Phi != NULL ) { 
    for (int v=0; v<Nband; v++) {
      for (int m=0; m<Npolar; m++ ) {
	for (int l=0; l<Nazim[m]; l++ ) {
	  for (int p=0; p<Npolar; p++ ) {
	    delete[] Phi[v][m][l][p];
	  } /* endfor */
	  delete[] Phi[v][m][l];
	} /* endfor */
	delete[] Phi[v][m];
      } /* endfor */
      delete[] Phi[v];
    } /* endfor */
    delete[] Phi; 
    Phi = NULL; 
  } /* endif */

}

/********************************************************
 * Object allocator and deallocator for SNBCK object.   *
 ********************************************************/
inline void Rte2D_State :: AllocateSNBCK()
{ 
  if ( SNBCKdata == NULL ) SNBCKdata = new SNBCK; 
}

inline void Rte2D_State :: DeallocateSNBCK()
{ 
  if ( SNBCKdata != NULL ) { delete SNBCKdata; SNBCKdata = NULL;}  
}

 /**************************************************************************
  ************************** OPERATOR OVERLOADS  ***************************
  **************************************************************************/

/********************************************************
 * Rte2D_State -- Binary arithmetic operators.          *
 ********************************************************/
// addition
inline Rte2D_State Rte2D_State :: operator +(const Rte2D_State &U) const {
  Rte2D_State Temp(*this);
  Temp += U;
  return Temp;
}

// subtraction
inline Rte2D_State Rte2D_State :: operator -(const Rte2D_State &U) const {
  Rte2D_State Temp(*this);
  Temp -= U;
  return Temp;
}

// scalar multiplication
inline Rte2D_State Rte2D_State :: operator *(const double &a) const {
  Rte2D_State Temp(*this);
  for ( int i=0; i<Ntot; i++ ) Temp.I[i] = a*I[i];
  return(Temp);
}

inline Rte2D_State operator *(const double &a, const Rte2D_State &U) {
  Rte2D_State Temp(U);
  for ( int i=0; i<U.Ntot; i++ ) Temp.I[i] = a*U.I[i];
  return(Temp);
}

// scalar division
inline Rte2D_State Rte2D_State :: operator /(const double &a) const {
  Rte2D_State Temp(*this);
  for ( int i=0; i<Ntot; i++ ) Temp.I[i]= I[i]/a;
  return(Temp);
}

// solution state division operator
inline Rte2D_State Rte2D_State :: operator /(const Rte2D_State &U) const {
  Rte2D_State Temp(*this);
  for ( int i=0; i<Ntot; i++ ) Temp.I[i]= I[i]/U.I[i];
  return(Temp);
}

// inner product
inline double Rte2D_State :: operator   *(const Rte2D_State &U) const {
  double sum=0.0;
  for ( int i=0; i<Ntot; i++ ) sum += I[i]*U.I[i];
  return sum;
}

// solution state product operator
inline Rte2D_State Rte2D_State :: operator ^(const Rte2D_State &U) const {
  Rte2D_State Temp(*this);
  for ( int i=0; i<Ntot; i++ ) Temp.I[i] = I[i]*U.I[i];
  return(Temp);
}

/********************************************************
 * Rte2D_State -- Assignment operator.                  *
 ********************************************************/
inline Rte2D_State& Rte2D_State :: operator =(const Rte2D_State &U) {
  if( this != &U) Copy(U);
  return (*this);
}

/********************************************************
 * Rte2D_State -- Shortcut operators.                   *
 ********************************************************/
inline Rte2D_State& Rte2D_State :: operator +=(const Rte2D_State &U) {
  for ( int i=0; i<Ntot; i++ ) I[i] += U.I[i];
  return (*this);
}

inline Rte2D_State& Rte2D_State :: operator -=(const Rte2D_State &U) {
  for ( int i=0; i<Ntot; i++ ) I[i] -= U.I[i];
  return (*this);
}

inline Rte2D_State& Rte2D_State::operator *=(const double &a) {
  for ( int i=0; i<Ntot; i++ ) I[i] *= a;
  return (*this);
}

inline Rte2D_State& Rte2D_State::operator /=(const double &a) {
  for ( int i=0; i<Ntot; i++ ) I[i] /= a;
  return (*this);
}

/********************************************************
 * Rte2D_State -- Unary arithmetic operators.           *
 ********************************************************/
inline Rte2D_State operator -(const Rte2D_State &U) {
  Rte2D_State Temp;
  for ( int i=0; i<U.Ntot; i++ ) Temp.I[i] = -U.I[i];
  return(Temp);
}

/********************************************************
 * Rte2D_State -- Relational operators.                 *
 ********************************************************/
inline int operator ==(const Rte2D_State &U1, const Rte2D_State &U2) {
  bool temp = true;
  for ( int i=0; i<U1.Ntot; i++ ) {
    if ( U1.I[i] != U2.I[i] ) { 
      temp = false; 
      break; 
    }
  }
  return temp;
}

inline int operator !=(const Rte2D_State &U1, const Rte2D_State &U2) {
  bool temp = false;
  for ( int i=0; i<U1.Ntot; i++ ) {
    if ( U1.I[i] != U2.I[i] ) { 
	  temp = true; 
	  break; 
    }
  }
  return temp;
}

inline int operator >(const Rte2D_State &U1, const double d) {
  bool temp = true;
  for ( int i=0; i<U1.Ntot; i++ ) {
    if ( U1.I[i] <= d ) { 
	  temp = false; 
	  break; 
    }
  }
  return temp;
}

inline int operator <(const Rte2D_State &U1, const double d) {
  bool temp = true;
  for ( int i=0; i<U1.Ntot; i++ ) {
    if ( U1.I[i] >= d ) { 
	  temp = false; 
	  break; 
    }
  }
  return temp;
}
 

/********************************************************
 * Rte2D_State -- Input/Output operators.               *
 ********************************************************/
inline ostream& operator << (ostream &out_file, const Rte2D_State &U) 
{
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.Ib[i];
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.kappa[i];
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.sigma[i];
  for( int i=0; i<U.Ntot; i++)  out_file<<" "<<U.I[i]; 
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream& operator >> (istream &in_file,  Rte2D_State &U) 
{
  in_file.setf(ios::skipws);
  for( int i=0; i<U.Nband; i++) in_file >> U.Ib[i];
  for( int i=0; i<U.Nband; i++) in_file >> U.kappa[i];
  for( int i=0; i<U.Nband; i++) in_file >> U.sigma[i];
  for( int i=0; i<U.Ntot; i++)  in_file >> U.I[i]; 
  in_file.unsetf(ios::skipws);
  return (in_file);
}

 /**************************************************************************
  ************************** FLUX/JAC FUNCTIONS  ***************************
  **************************************************************************/

/*************************************************************
 * Rte2D_State::Fn -- Flux in n-direction (Time March).      *
 *************************************************************/
inline Rte2D_State Rte2D_State :: Fn( const Vector2D &norm_dir ) const 
{

  double cos_angle, sin_angle, cosine;
  Rte2D_State Temp;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;


  for (int v=0; v<Nband; v++ ) 
    for (int m=0; m<Npolar; m++) 
      for (int l=0; l<Nazim[m]; l++) {
	cosine = mu[m][l]*cos_angle + eta[m][l]*sin_angle;
	Temp.In(v,m,l) = cosine*In(v,m,l);
      }
  
  return (Temp);
}

/*************************************************************
 * Rte2D_State::Fn -- Flux in n-direction (Space March).     *
 *************************************************************/
inline double Rte2D_State :: Fn( const Vector2D &norm_dir, const int v, 
				 const int m, const int l ) const 
{

  double cos_angle, sin_angle, cosine;
  Rte2D_State Temp;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;
  cosine = mu[m][l]*cos_angle + eta[m][l]*sin_angle;

  return (cosine*In(v,m,l));
}


/*************************************************************
 * Rte2D_State::dFdU -- Flux Jacobian in n-direction.        *
 *                      Time-march version.                  *
 *************************************************************/
inline void dFndU(DenseMatrix &dFdU, const Rte2D_State &U, const Vector2D &norm_dir) {

  double cos_angle, sin_angle, cosine;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  int index = 0;
  for (int v=0; v<U.Nband; v++ ) 
    for (int m=0; m<U.Npolar; m++) 
      for (int l=0; l<U.Nazim[m]; l++) {
	cosine = U.mu[m][l]*cos_angle + U.eta[m][l]*sin_angle;

// 	cout.precision(4);
// 	cout << setw(4) << m+1 
// 	     << setw(4) << l+1
// 	     << setw(10) << U.mu[m][l]
// 	     << setw(10) << U.eta[m][l]
// 	     << setw(10) << U.mu[m][l]*cos_angle + U.eta[m][l]*sin_angle
// 	     << endl;
	if ( cosine < ZERO )
	  dFdU(index,index) +=  cosine;
	index++;
      }
}


/*************************************************************
 * Rte2D_State::dFdU -- Flux Jacobian in n-direction.        *
 *                      Space-march version.                 *
 *************************************************************/
inline double Rte2D_State :: dFndU(const Vector2D &norm_dir, 
				   const int m, const int l)  const {

  double cos_angle, sin_angle, cosine;

  /* Determine the direction cosine's for the frame rotation. */
  cos_angle = norm_dir.x; 
  sin_angle = norm_dir.y;

  cosine = mu[m][l]*cos_angle + eta[m][l]*sin_angle;

  return cosine;
}




/*************************************************************
 * Rte2D_State::s -- Regular source term.                    *
 *                   Time-march version.                     *
 *************************************************************/
inline Rte2D_State Rte2D_State :: S(void) const 
{
  // declares
  double temp;
  double beta;
  Rte2D_State Temp;
  Temp.Zero();

  // loop along bands
  for (int v=0; v<Nband; v++ ) {

    // the extinction coefficient
    beta = kappa[v] + sigma[v];
    
    // loop over directions
    for (int m=0; m<Npolar; m++) 
      for (int l=0; l<Nazim[m]; l++) {

	// add blackbody emission
	Temp.In(v,m,l)  = kappa[v] * Ib[v];

	// subtract absorbsion and out-scattering
	Temp.In(v,m,l) -= beta * In(v,m,l);

	// add in-scattering
	if (sigma[v]>TOLER) {
	  temp = ZERO;
	  for (int p=0; p<Npolar; p++) 
	    for (int q=0; q<Nazim[p]; q++) {
	      temp += In(v,p,q) * Phi[v][p][q][m][l] * omega[p][q];
	    } /* endfor - in-dirs*/
	  temp *= sigma[v] / (FOUR * PI);
	  Temp.In(v,m,l) += temp;
	} /* endif - in-scat */
	  
	// multiply
	Temp.In(v,m,l) *= omega[m][l];

      } /* endfor - dirs*/

  } /* endfor - bands */

  return (Temp);
}


/*************************************************************
 * Rte2D_State::s -- Regular source term.                    *
 *                   Space-march version.                    *
 *************************************************************/
inline double Rte2D_State :: S(const int v, const int m, const int l) const 
{
  // declares
  double temp2;
  double beta;
  double temp = ZERO;


  // the extinction coefficient
  beta = kappa[v] + sigma[v];
    
  // add blackbody emission
  temp = kappa[v] * Ib[v];

  // add in-scattering
  if (sigma[v]>TOLER) {
    temp2 = ZERO;
    for (int p=0; p<Npolar; p++) 
      for (int q=0; q<Nazim[p]; q++) {
	// skip forward scattering, it is acoundted for in dSdU
	if(p==m && l==q) continue;
	temp2 += In(v,p,q) * Phi[v][p][q][m][l] * omega[p][q];
      } /* endfor - in-dirs*/
    temp2 *= sigma[v] / (FOUR * PI);
    temp += temp2;
  } /* endif - in-scat */
	  
  // multiply
  temp *= omega[m][l];


  return (temp);
}


/*************************************************************
 * Rte2D_State:: dSdU -- Regular source term jacobian.       *
 *                       Time-march version.                 *
 *************************************************************/
inline void Rte2D_State :: dSdU(DenseMatrix &dSdU) const
{ 
  // declares
  double temp;
  double beta;

  // loop along bands
  for (int v=0; v<Nband; v++ ) {

    // the extinction coefficient
    beta = kappa[v] + sigma[v];
    
    // loop over directions
    for (int m=0; m<Npolar; m++) 
      for (int l=0; l<Nazim[m]; l++) {


	// absorbsion and out-scattering
	dSdU(Index[v][m][l],Index[v][m][l]) -= omega[m][l] * beta;

	// in-scattering
	if (sigma[v]>TOLER) {

	  temp = sigma[v] * omega[m][l] / (FOUR * PI);
	  for (int p=0; p<Npolar; p++) 
	    for (int q=0; q<Nazim[p]; q++){
	      dSdU(Index[v][m][l],Index[v][p][q]) += Phi[v][p][q][m][l] * omega[p][q] * temp;
	    }

	} /* endif - in-scat */
	  
      } /* endfor - dirs*/

  } /* endfor - bands */

}

/*************************************************************
 * Rte2D_State:: dSdU -- Regular source term jacobian.       *
 *                       Space-march version.                *
 *************************************************************/
inline double Rte2D_State :: dSdU(const int v, const int m, const int l) const
{ 
  // declares
  double temp = ZERO;
  double beta;

  // the extinction coefficient
  beta = kappa[v] + sigma[v];
    
  // absorbsion and out-scattering
  temp = omega[m][l] * beta;

  // in-scattering
  if (sigma[v]>TOLER) {
    
    // add forward scattering
    temp -= sigma[v] * omega[m][l] / (FOUR * PI) 
      * Phi[v][m][l][m][l] * omega[m][l];
    
  } /* endif - in-scat */
  
  return temp;

}



/*************************************************************
 * Rte2D_State::s_axi -- Axisymmetric source term.           *
 *                       Time-march version.                 *
 *************************************************************/
inline Rte2D_State Rte2D_State :: Sa(const Rte2D_State &dUdpsi,
				     const Rte2D_State &phi_psi,
				     const int Axisymmetric) const 
{
  int l_m1, l_p1;
  double It, Ib, Dt, Db;
  Rte2D_State Temp;  Temp.Zero();

  // assume constant angular spacing
  // in azimuthal angle
  
  // loop along bands and directions

  for (int v=0; v<Nband; v++ ) 
    for (int m=0; m<Npolar; m++) 
      for (int l=0; l<Nazim[m]; l++) {
	
	// use symmetry to get l-1 and l+1
	if (l==0) l_m1 = l;
	else l_m1 = l-1;
	
	if (l==Nazim[m]-1) l_p1 = l;
	else l_p1 = l+1;
	 
	// direction cosines
	if (Axisymmetric == AXISYMMETRIC_Y) {
	  Dt = -sin(delta_psi[m][l]/TWO)*eta[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
	  Db = -sin(delta_psi[m][l]/TWO)*eta[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];    
	} else if (Axisymmetric == AXISYMMETRIC_X) {
	  Dt = -sin(delta_psi[m][l]/TWO)*mu[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
	  Db = -sin(delta_psi[m][l]/TWO)*mu[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];
	} /* endif */
	
	// upwind
	if (Dt<ZERO) {
	  It = In(v,m,l_m1) +
	    HALF*phi_psi.In(v,m,l_m1)*dUdpsi.In(v,m,l_m1)*delta_psi[m][l_m1];
	} else {
	  It = In(v,m,l) - 
	    HALF*phi_psi.In(v,m,l)*dUdpsi.In(v,m,l)*delta_psi[m][l];
	}
	if (Db>ZERO) {
	  Ib = In(v,m,l) + 
	    HALF*phi_psi.In(v,m,l)*dUdpsi.In(v,m,l)*delta_psi[m][l];
	} else {
	  Ib = In(v,m,l_p1) - 
	    HALF*phi_psi.In(v,m,l_p1)*dUdpsi.In(v,m,l_p1)*delta_psi[m][l_p1];
	}

	// compute source term
	Temp.In(v,m,l) -= ( Dt*It + Db*Ib )/delta_psi[m][l];
		
	
      } /* endfor - dirs*/

  return (Temp);
}

/*************************************************************
 * Rte2D_State::s_axi -- Axisymmetric source term.           *
 *                       Space-march FVM version.            *
 *************************************************************/
inline double Rte2D_State :: Sa_FVM(const int v, const int m, const int l,
				    const int Axisymmetric) const 
{
  int l_m1, l_p1;
  double It, Ib, Dt, Db;
  double Temp;

  // assume constant angular spacing
  // in azimuthal angle
  
  // use symmetry to get l-1 and l+1
  if (l==0) l_m1 = l;
  else l_m1 = l-1;
  
  if (l==Nazim[m]-1) l_p1 = l;
  else l_p1 = l+1;
  
  // direction cosines
  if (Axisymmetric == AXISYMMETRIC_Y) {
    Dt = -sin(delta_psi[m][l]/TWO)*eta[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
    Db = -sin(delta_psi[m][l]/TWO)*eta[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];    
  } else if (Axisymmetric == AXISYMMETRIC_X) {
    Dt = -sin(delta_psi[m][l]/TWO)*mu[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
    Db = -sin(delta_psi[m][l]/TWO)*mu[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];
  } /* endif */
  
  // upwind
  if (Dt<ZERO) {
    It = In(v,m,l_m1);
  } else {
    It = ZERO; // It = In(v,m,l); -> goes into dSadU term
  }
  if (Db>ZERO) {
    Ib = ZERO; // Ib = In(v,m,l); -> goes into dSadU term
  } else {
    Ib = In(v,m,l_p1);
  }
  
  // compute source term
  Temp = -( Dt*It + Db*Ib )/delta_psi[m][l];
  	
	

  return (Temp);
}


/*************************************************************
 * Rte2D_State::s_axi -- Axisymmetric source term.           *
 *                       Space-march DOM version.            *
 *************************************************************/
inline double Rte2D_State :: Sa_DOM(const int v, const int m, const int l) const 
{
  int l_m1, l_p1;
  double Temp;

  // use symmetry to get l-1 and l+1
  if (l==0) l_m1 = l;
  else l_m1 = l-1;
  
  if (l==Nazim[m]-1) l_p1 = l;
  else l_p1 = l+1;
  
  // compute source term
  Temp = -alpha[m][l]*I_half[m][l];
  	
	

  return (Temp);
}


/*************************************************************
 * Rte2D_State::dSadU -- Axisymmetric source term jacobian.  *
 *                       Time-march version.                 *
 *************************************************************/
inline void Rte2D_State :: dSadU(DenseMatrix &dSadU,
				 const double Sp,
				 const int Axisymmetric) const 
{ 

  int l_m1, l_p1;
  double Dt, Db;

  // loop along bands
  for (int v=0; v<Nband; v++ ) {

   
    // loop over directions
    for (int m=0; m<Npolar; m++) 
      for (int l=0; l<Nazim[m]; l++) {


	// use symmetry to get l-1 and l+1
	if (l==0) l_m1 = l;
	else l_m1 = l-1;
	
	if (l==Nazim[m]-1) l_p1 = l;
	else l_p1 = l+1;
	 
	// direction cosines
	if (Axisymmetric == AXISYMMETRIC_Y) {
	  Dt = -sin(delta_psi[m][l]/TWO)*eta[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
	  Db = -sin(delta_psi[m][l]/TWO)*eta[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];    
	} else if (Axisymmetric == AXISYMMETRIC_X) {
	  Dt = -sin(delta_psi[m][l]/TWO)*mu[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
	  Db = -sin(delta_psi[m][l]/TWO)*mu[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];
	} /* endif */
	
	// upwind
	if (Dt<ZERO) {
	  // It = I[ Index[v][m][l_m1] ];
	  dSadU(Index[v][m][l],Index[v][m][l_m1]) -= Dt/(delta_psi[m][l]*Sp);
	} else {
	  // It = I[ Index[v][m][l] ];
	  dSadU(Index[v][m][l],Index[v][m][l]) -= Dt/(delta_psi[m][l]*Sp);
	}
	if (Db>ZERO) {
	  // Ib = I[ Index[v][m][l] ];
	  dSadU(Index[v][m][l],Index[v][m][l]) -= Db/(delta_psi[m][l]*Sp);
	} else {
	  // Ib = I[ Index[v][m][l_p1] ];
	  dSadU(Index[v][m][l],Index[v][m][l_p1]) -= Db/(delta_psi[m][l]*Sp);
	}
	

      } /* endfor - dirs*/

  } /* endfor - bands */


}


/*************************************************************
 * Rte2D_State::dSadU -- Axisymmetric source term jacobian.  *
 *                       Space-march FVM version.            *
 *************************************************************/
inline double Rte2D_State :: dSadU_FVM(const int v, const int m, const int l, 
				       const int Axisymmetric) const 
{ 

  int l_m1, l_p1;
  double Dt, Db;
  double Temp;


  // use symmetry to get l-1 and l+1
  if (l==0) l_m1 = l;
  else l_m1 = l-1;
  
  if (l==Nazim[m]-1) l_p1 = l;
  else l_p1 = l+1;
  
  // direction cosines
  if (Axisymmetric == AXISYMMETRIC_Y) {
    Dt = -sin(delta_psi[m][l]/TWO)*eta[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
    Db = -sin(delta_psi[m][l]/TWO)*eta[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];
  } else if (Axisymmetric == AXISYMMETRIC_X) {
    Dt = -sin(delta_psi[m][l]/TWO)*mu[m][l] + cos(delta_psi[m][l]/TWO)*xi[m][l];
    Db = -sin(delta_psi[m][l]/TWO)*mu[m][l] - cos(delta_psi[m][l]/TWO)*xi[m][l];
  } /* endif */
  
  // upwind
  if (Dt<ZERO) {
    // It = I[ Index[v][m][l_m1] ];
    Temp = ZERO;
  } else {
    // It = I[ Index[v][m][l] ];
    Temp = Dt/(delta_psi[m][l]);
  }
  if (Db>ZERO) {
    // Ib = I[ Index[v][m][l] ];
    Temp += Db/(delta_psi[m][l]);
  } else {
    // Ib = I[ Index[v][m][l_p1] ];
    Temp += ZERO;
  }
  

  return Temp;

}


/*************************************************************
 * Rte2D_State::dSadU -- Axisymmetric source term jacobian.  *
 *                       Space-march DOM version.            *
 *************************************************************/
inline double Rte2D_State :: dSadU_DOM(const int v, const int m, const int l) const 
{ 

  int l_m1, l_p1;
  double Temp;


  // use symmetry to get l-1 and l+1
  if (l==0) l_m1 = l;
  else l_m1 = l-1;
  
  if (l==Nazim[m]-1) l_p1 = l;
  else l_p1 = l+1;
  
  
  // upwind
  Temp = - alpha[m][l+1];

  return Temp;

}



 /**************************************************************************
  ************************** EXTERNAL FUNCTIONS  ***************************
  **************************************************************************/

/*******************************************************
 * External Flux Function Functions                    *
 *******************************************************/  
extern Rte2D_State Riemann_n(const Rte2D_State &Ul,
			     const Rte2D_State &Ur,
			     const Vector2D &norm_dir );

extern Rte2D_State Flux_n(const Rte2D_State &Ul,
			  const Rte2D_State &Ur,
			  const Vector2D &norm_dir);

/********************************************************
 * External Boundary Conditions Functions               *
 ********************************************************/
extern Rte2D_State Reflect(const Rte2D_State &U,
			   const Vector2D &norm_dir);
extern void Reflect_Space_March(Rte2D_State &U,
				const Vector2D &norm_dir);

extern Rte2D_State Gray_Wall(const Rte2D_State &U, 
		      const Vector2D &norm_dir, 
		      const double &wall_temperature,
		      const double &wall_emissivity );
extern void Gray_Wall_Space_March(Rte2D_State &Uwall, 
		      const Vector2D &norm_dir, 
		      const double &wall_temperature,
		      const double &wall_emissivity );

/********************************************************
 * Miscellaneous Scattering Functions                   *
 ********************************************************/
extern double Legendre( const double x, const int n);
extern double* PhaseFunc( const int type, int &n);

#endif //end _RTE2D_STATE_INCLUDED 
