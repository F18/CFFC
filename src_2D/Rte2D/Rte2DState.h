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

// Required C++ libraries
#include <iostream>
using namespace std;

// Required CFDkit+caboodle header files

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

#ifndef _VECTOR3D_INCLUDED
#include "../Math/Vector3D.h"
#endif //_VECTOR3D_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _SIMPSON_INCLUDED
#include "../Math/Simpson.h"
#endif // _SIMPSON_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED


// Constants
/********************************************************
 * Gas Type                                             *
 ********************************************************/
enum Gas_Models { RTE2D_GAS_TYPE_GRAY };

/********************************************************
 * RTE discretization Type (DOM or FVM)                 *
 ********************************************************/
enum RTE_Discretization { RTE2D_SOLVER_FVM,
			  RTE2D_SOLVER_DOM };

/********************************************************
 * Scatter Phase Functions                              *
 ********************************************************/
enum Scatter_Models { SCATTER_ISO,  // isotropic scattering
		      SCATTER_F1,   // Kim and Lee (1988)   
		      SCATTER_F2,   // Kim and Lee (1988)
		      SCATTER_F3,   // Kim and Lee (1988) 
		      SCATTER_B1,   // Kim and Lee (1988)
		      SCATTER_B2 }; // Kim and Lee (1988)

/********************************************************
 * DOM quadrature Type                                  *
 ********************************************************/
enum DOM_Qauad { DOM_S2,   // S2 NONSYMMETRIC from LATHROP and CARLSON
		 DOM_S4,   // S4 SYMMETRIC from LATHROP and CARLSON
		 DOM_S6,   // S6 SYMMETRIC from LATHROP and CARLSON
		 DOM_S8,   // S8 SYMMETRIC from LATHROP and CARLSON
		 DOM_S12,  // S12 SYMMETRIC from LATHROP and CARLSON
		 DOM_T3 }; // T3 SYMMETRIC from Truelove


//-------------- STRUCTS REQUIRED FOR INTEGRATION -------------------//

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

//-------------- END STRUCTS REQUIRED FOR INTEGRATION ----------------//




/**********************************************************************
 * Class: Rte2D_State                                                    *
 *                                                                    *
 * Description: The radiation state class contains the properties     *
 *              of a gray, absorbing, emitting, anisotropically       *
 *              scattering medium.  Additionally, it contains the     *
 *              state (intensity) of a beam of light passing          *
 *              through the medium.                                   *
 *                                                                    *
 **********************************************************************/

class Rte2D_State {

  
  //---------------------------- Objects ----------------------------//
  
 public:

  double* kappa;              // absorbsion coefficient [m^-1]
  double* sigma;              // scattering coefficient [m^-1]
  double* Ib;                 // blackbody intentsity [W/m^2]
  double* I;                  // directional intentsity [W/m^2]
  static int  Npolar;         // number of polar directions
  static int* Nazim;          // number of azimuthal directions
  static int  Nband;          // the total number of frequency bands
  static int  Ntot;           // total number of intensities
  static int*** Index;        // indexing array
  static int NUM_VAR_RTE2D;   // total number of Rte2D variables
  static double** mu;         // (x,r)-direction cosine
  static double** eta;        // (y,azimuthal)-direction cosine
  static double** xi;         // (z)-direction cosine
  static double** omega;      // discrete control angle element sizes
  static double* theta;       // polar angle grid points
  static double** psi;        // azimuthal angle grid points
  static double* delta_theta; // polar angle grid points
  static double** delta_psi;  // azimuthal angle grid points
  static double***** Phi;     // scattering phase function
  static int Symmetry_Factor; // solid angle range symmetry factor
  static int RTE_Type;        // flag for DOM or FVM
  double** I_half;             // storage for intensity in special directions
                              // See Carlson and Lathrop (1968)

 private:

  double** alpha;             // angular redistribution coefficient (for DOM)

  //--------------------------- Functions ---------------------------//

 public:

  // constructors
  Rte2D_State();
  Rte2D_State( const Rte2D_State &U );
  
  // destructors
  ~Rte2D_State();

  // copy function
  void Copy( const Rte2D_State &U );
  void Copy_NonSol( const Rte2D_State &U );

  // memory allocation / deallocation for the I array
  void Allocate();
  void Deallocate();

  // memory allocation / deallocation for the directional arrays
  static void AllocateDirs( );
  static void DeallocateDirs();

  // memory allocation / deallocation for the cosine arrays
  static void AllocateCosines( );
  static void DeallocateCosines();

  // Zero the state
  void Zero();
  void zero_non_sol();

  // check for negative intensities
  int NegIntensity();

  // initialize
  void SetIntensity(double val);
  void SetAbsorbsion(double val);
  void SetScattering(double val);
  void SetBlackbody(double val);

  // access intensity using 2D indexes 
  double& In( const int v, const int m, const int l ) const;
  double& In( const int v, const int m, const int l );

  // compute extinction coefficient
  double beta(const int v) const;

  // set the number of directions and compute the cosines
  static void SetDirs(const int NumPolarDirs, const int NumAzimDirs, 
		      const int Quad_Type, const int Axisymmetric);
  static void SetDirsFVM(const int NumPolarDirs, const int NumAzimDirs, 
			 const int Quad_Type, const int Axisymmetric);
  static void SetDirsDOM(const int NumPolarDirs, const int NumAzimDirs, 
			 const int Quad_Type, const int Axisymmetric);

  // setup the phase function
  static void SetupPhase( const int type );
  static void SetupPhaseFVM( const int type );
  static void SetupPhaseDOM( const int type );

  // set the number of directions and compute the cosines
  static void SetGas();

  // Flux.
  Rte2D_State Fn( const Vector2D &norm_dir ) const;
  double Fn( const Vector2D &norm_dir, const int v, 
	     const int m, const int l ) const;
  // Rte2D_State Fx(void) const;
  // Rte2D_State Fy(void) const;

  // Flux Jacobian X 
  friend void dFdU(DenseMatrix &dFdU, const Rte2D_State &U);
  friend void dFdU_n(DenseMatrix &dFdU, const Rte2D_State &U, 
		     const Vector2D &norm_dir);
  double dFdU_n( const Vector2D &norm_dir, const int m, const int l ) const;
  Rte2D_State dFdU_n(const Vector2D &norm_dir) const;


  // Regular source term.
  Rte2D_State S(void) const;
  double S(const int v, const int m, const int l) const;

  // source term jacobian
  void dSdU(DenseMatrix &dSdU) const;
  double dSdU(const int v, const int m, const int l ) const;

  // Axisymmetric source term. 
  Rte2D_State Sa( const Rte2D_State &dUdpsi,
		  const Rte2D_State &phi_psi ) const;
  double Sa_FVM(const int v, const int m, const int l) const;
  double Sa_DOM(const int v, const int m, const int l) const;

  // Axisymmetric source term jacobian
  void dSadU(DenseMatrix &dSadU, const double Sp) const;
  double dSadU_FVM(const int v, const int m, const int l) const;
  double dSadU_DOM(const int v, const int m, const int l) const;

  // compute DOM angular redistribution coefficients
  void SetupART_DOM( const Vector2D nfaceE, const double AfaceE,
		     const Vector2D nfaceW, const double AfaceW,
		     const Vector2D nfaceN, const double AfaceN,
		     const Vector2D nfaceS, const double AfaceS );


  // Compute directional integrated intensity [W/m^2]
  double G();        

  // heat flux vector [W/m^2]
  Vector2D q();


  //--------------------------- Operators ---------------------------//

  // Index operator 
  double &operator[](int index);
  const double &operator[](int index) const;

  // Binary arithmetic operators.
  Rte2D_State operator +(const Rte2D_State &U) const;
  Rte2D_State operator -(const Rte2D_State &U) const;
  Rte2D_State operator *(const double &a) const;
  double operator   *(const Rte2D_State &U) const;
  friend Rte2D_State operator *(const double &a, const Rte2D_State &U);
  Rte2D_State operator /(const double &a) const;
  Rte2D_State operator /(const Rte2D_State &U) const;
  Rte2D_State operator ^(const Rte2D_State &U) const;

  // Assignment Operator. 
  Rte2D_State& operator =(const Rte2D_State &U); 

  // Shortcut arithmetic operators.
  Rte2D_State& operator +=(const Rte2D_State &U);
  Rte2D_State& operator -=(const Rte2D_State &U);
      
  // Unary arithmetic operators. 
  friend Rte2D_State operator -(const Rte2D_State &U);
  
  // Relational operators.
  friend int operator ==(const Rte2D_State &U1, const Rte2D_State &U2);
  friend int operator !=(const Rte2D_State &U1, const Rte2D_State &U2);
  friend int operator >(const Rte2D_State &U1, const double d);
  friend int operator <(const Rte2D_State &U1, const double d);

  // Input-output operators.
  friend ostream& operator << (ostream &out_file, const Rte2D_State &U);
  friend istream& operator >> (istream &in_file,  Rte2D_State &U);



};


/********************************************************
 *                       EXTERNAL FUNCTIONS             *
 ********************************************************/

/********************************************************
 * Blackbody intensity                                  *
 ********************************************************/
extern double Ib(const double T);
extern double Ib(const double T, const double lambda1, 
		 const double lambda2);
extern double Planck(const double lambdaT);


/*******************************************************
 * External Flux Function Functions                    *
 *******************************************************/  
// extern Rte2D_State Riemann(const Rte2D_State &Ul,
// 			   const Rte2D_State &Ur);

// extern Rte2D_State Riemann_x(const Rte2D_State &Ul,
// 			     const Rte2D_State &Ur);

// extern Rte2D_State Riemann_y(const Rte2D_State &Ul,
// 			     const Rte2D_State &Ur);
extern Rte2D_State Riemann_n(const Rte2D_State &Ul,
			     const Rte2D_State &Ur,
			     const Vector2D &norm_dir );

// extern Rte2D_State Flux(const Rte2D_State &Ul,
// 			const Rte2D_State &Ur);
// extern Rte2D_State Flux_x(const Rte2D_State &Ul,
// 			  const Rte2D_State &Ur);
// extern Rte2D_State Flux_y(const Rte2D_State &Ul,
// 			  const Rte2D_State &Ur);
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
extern double func(int ndim, double *x, void *params);

/********************************************************
 * Multigrid restriction and prolongation functions     *
 ********************************************************/
extern void Restrict_NonSol( Rte2D_State &Uc, const double Ac,
			     const Rte2D_State &Uf_SW, const double Af_SW,
			     const Rte2D_State &Uf_SE, const double Af_SE,
			     const Rte2D_State &Uf_NW, const double Af_NW,
			     const Rte2D_State &Uf_NE, const double Af_NE );
extern void Restrict_NonSol_Boundary_Ref_States( 
			     Rte2D_State &Uc, 
			     const Rte2D_State &Uf_l, const double Af_l,
			     const Rte2D_State &Uf_r, const double Af_r );


extern int Prolong_NonSol( const Rte2D_State &Uc1, const Vector2D Xc1,
			   const Rte2D_State &Uc2, const Vector2D Xc2,
			   const Rte2D_State &Uc3, const Vector2D Xc3,
			   const Rte2D_State &Uc4, const Vector2D Xc4,
			   const Vector2D XcP, Rte2D_State &Uf );


/********************************************************
 * Exact solution functions                             *
 ********************************************************/
extern void CylindricalEnclosure( const double gas_temp,
				  const double c,
				  const double tau,
				  const double rpos,
				  const double zpos,
				  double &G, 
				  double &qr, 
				  double &qz );

extern void RectangularEnclosure( const double gas_temp,
				  const double kappa,
				  const double left,
				  const double right,
				  const double bot,
				  const double top,
				  const double xpos,
				  const double ypos,
				  double &G, 
				  double &qx, 
				  double &qy );



#endif //end _RTE2D_STATE_INCLUDED 
