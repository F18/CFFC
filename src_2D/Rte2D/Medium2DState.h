/**********************************************************************
 **********************************************************************
 **                                                                  **
 ** File: Medium2DState.h                                            **
 **                                                                  **
 ** Description: The medium state class contains the properties      **
 **              of a gray, absorbing, emitting, anisotropically     **
 **              scattering medium.   This file defines the state    ** 
 **              class.                                              **
 **                                                                  **
 ** Author: Marc "T-Bone" Charest                                    **
 **                                                                  **
 ** Revision:  Date        Initials   Change                         **
 **            04/03/2007  MRC        Original creation              **
 **                                                                  **
 **********************************************************************
 **********************************************************************/
#ifndef _Medium2D_STATE_INCLUDED
#define _Medium2D_STATE_INCLUDED 

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
 * Required CFFC header files                           *
 ********************************************************/

#include "../Math/Math.h"
#include "../Math/Vector2D.h"
#include "../CFD/CFD.h"
#include "SNBCK.h"
#include "Vector2D_Function.h"

/********************************************************
 * Necessary Medium2D Specific Constants                *
 ********************************************************/

// Absorbsion model type 
enum Gas_Models { MEDIUM2D_ABSORB_GRAY,
                  MEDIUM2D_ABSORB_SNBCK };

// Absorbsion model type 
enum Field_Type { MEDIUM2D_FIELD_ANALYTIC,
                  MEDIUM2D_FIELD_DISCRETE };

/***********************************************************************/
/*!
 * Class: Medium_State
 *
 * @brief Gas state class definition for a absorbing, scattering 
 *        participating medium.
 *
 * Gas state class definition for a absorbing, scattering participating
 * medium.  This is a container for prescirbed field data required by
 * the Rte2D class.
 *
 * \verbatim
 * Member functions
 *     kappa    -- Return array of spectral absorbsion coefficient
 *     sigma    -- Return array of spectral scattering coefficient
 *     Ib       -- Return array of spectral blackbody intensity.
 *     beta     -- Return extinction coefficient.
 *     Nbands   -- Return the number of spectral frequency bands.
 * \endverbatim
 */
/***********************************************************************/
class Medium2D_State {

 private:
  
  //@{ @name Field data:
  static Vector2D_Function<Medium2D_State>* Field;  //!< container class for analytic function
  //@}


 public:

  //@{ @name medium state:
  double* kappa;        //!< absorbsion coefficient [m^-1]
  double* sigma;        //!< scattering coefficient [m^-1]
  double* Ib;           //!< blackbody intentsity [W/m^2 or W/(m^2 cm)]
  //@}

  //@{ @name Public objects:
  static int Nband;            //!< the total number of frequency bands (and quadrature points for SNBCK)
  static int NUM_VAR_MEDIUM2D; //!< the total number of variables in the state
  static SNBCK* SNBCKdata;     //!< statistical narrow band model 
  static double Absorb_Type;   //!< flag for absorption model
  //@}


  //@{ @name Creation, copy, and assignment constructors.
  //! Creation constructor.
  Medium2D_State() : kappa(NULL), sigma(NULL), Ib(NULL)
    { Allocate(); }

  //! Copy constructor.
  Medium2D_State( const Medium2D_State &U ) : kappa(NULL), sigma(NULL), Ib(NULL)
    { Allocate(); if( this != &U) Copy(U); }
  
  //! Destructor.
  ~Medium2D_State() { Deallocate(); }
  //@}

  //@{ @name Useful operators.
  //! Copy non solution state operator.
  void Copy( const Medium2D_State &U );

  //! Zero operator.
  void Zero();  

  //! Initialer
  void SetInitialValues( const double &Pressure,
			 const double &Temperature,
			 const double &xco,
			 const double &xh2o,
			 const double &xco2,
			 const double &xo2,
			 const double &fsoot,
			 const double &AbsorptionCoef,
			 const double &ScatteringCoef);
  //@}

  //@{ @name Allocators and deallocators
  //! memory allocation / deallocation for the private arrays
  void Allocate();
  void Deallocate();
  //! memory allocation / deallocation for the SNBCK data object
  static void AllocateSNBCK();
  static void DeallocateSNBCK();
  //! memory allocation / deallocation for the field object
  static void DeallocateField();
  //! deallocate all static variables
  static void DeallocateAllStatic() { DeallocateSNBCK(); DeallocateField(); }
  //@}

  //@{ @name State functions.
  //! Return extinction coefficient
  double beta(const int &v) const { return (kappa[v] + sigma[v]); }

  //! Compute analytic medium state at location
  void SetState(const Vector2D &r);
  static Medium2D_State GetState(const Vector2D &r);

  //! Compute a new state dependant upon gas state
  void ComputeNewState( const double &Pressure,
			const double &Temperature,
			const double &xco,
			const double &xh2o,
			const double &xco2,
			const double &xo2,
			const double &fsoot );
  //@}


  //@{ @name Static functions
  //! Setup function
  static void SetupStatic( const int &i_Absorb_Type, 
			   const SNBCK_Input_Parameters &SNBCK_IP,
			   const char* PATH);

  //! Set all all fields to the same function
  static void SetConstantField(const Medium2D_State &M);
  //@}

  //@{ @name Assignment operator.
  Medium2D_State& operator =(const Medium2D_State &U) {
    if( this != &U) Copy(U);
    return (*this);
  }
  //@}


  //@{ @name Index operator.
  double& operator[](int index);
  const double& operator[](int index) const;
  //@}

  //@{ @name Arithmatic operators.
  //! Binary arimatic operators
  Medium2D_State operator -(const Medium2D_State &U) const;
  Medium2D_State operator +(const Medium2D_State &U) const;
  Medium2D_State operator *(const double &a) const;
  friend Medium2D_State operator *(const double &a, const Medium2D_State &U);
  Medium2D_State operator /(const double &a) const;
 //!Shortcut operators
  Medium2D_State& operator -=(const Medium2D_State &U);
  Medium2D_State& operator +=(const Medium2D_State &U);
  Medium2D_State& operator *=(const double &a);
  Medium2D_State& operator /=(const double &a);
  //@}


  //@{ @name Input-output operators.
  friend ostream& operator << (ostream &out_file, const Medium2D_State &U);
  friend istream& operator >> (istream &in_file,  Medium2D_State &U);
 //@}



};

/********************************************************
 * Copy function.                                       *
 ********************************************************/
inline void Medium2D_State :: Copy( const Medium2D_State &U ) {
  for ( int i=0; i<Nband; i++ )   { 
    Ib[i] = U.Ib[i];
    kappa[i] = U.kappa[i];
    sigma[i] = U.sigma[i];
  }
}

/********************************************************
 * Zero operator.                                       *
 ********************************************************/
inline void Medium2D_State :: Zero() {
  for(int i=0; i<Nband; i++) {
    kappa[i] = ZERO;  
    sigma[i] = ZERO;
    Ib[i]    = ZERO;  
  }
}  

/********************************************************
 * Array allocator and deallocator for dynamic arrays.  *
 ********************************************************/
inline void Medium2D_State :: Allocate()
{
  // deallocate first
  Deallocate();

  // create the arrays
  if (Nband>0) {
    Ib = new double[Nband];
    kappa = new double[Nband];
    sigma = new double[Nband];
  }
}

inline void Medium2D_State :: Deallocate()
{
  if ( kappa != NULL ) { delete[] kappa;   kappa = NULL; }
  if ( sigma != NULL ) { delete[] sigma;   sigma = NULL; }
  if (    Ib != NULL ) { delete[] Ib;         Ib = NULL; }
}


inline void Medium2D_State :: AllocateSNBCK()
{ 
  DeallocateSNBCK();
  SNBCKdata = new SNBCK; 
}

inline void Medium2D_State :: DeallocateSNBCK()
{ if ( SNBCKdata != NULL ) { delete SNBCKdata; SNBCKdata = NULL;} }


inline void Medium2D_State :: DeallocateField()
{ if ( Field != NULL ) { delete Field; Field = NULL;}  }



/********************************************************
 * Compute medium state at location.                    *
 ********************************************************/
inline void Medium2D_State :: SetState(const Vector2D &r) {
  *this = (*Field)(r);
}


inline Medium2D_State Medium2D_State :: GetState(const Vector2D &r) {
  return (*Field)(r);
}


/********************************************************
 * Index operator.                                      *
 ********************************************************/
// WARNING - this should be used sparingly, it is not very efficient
inline double& Medium2D_State :: operator[](int index){
  if (index>0 && index<=Nband) return kappa[index-1];
  else if (index>Nband && index<=2*Nband) return sigma[index-Nband-1];
  else if (index>2*Nband && index<=3*Nband) return Ib[index-2*Nband-1];
  else {
    cerr << "Medium2D_State::operator[] - out of bounds.\n";
    exit(-1);
  }
}

inline const double& Medium2D_State :: operator[](int index) const {
  if (index<=Nband) return kappa[index-1];
  else if (index>Nband && index<=2*Nband) return sigma[index-Nband-1];
  else if (index>2*Nband && index<=3*Nband) return Ib[index-2*Nband-1];
  else {
    cerr << "Medium2D_State::operator[] - out of bounds.\n";
    exit(-1);
  }
}


/********************************************************
 * Binary arithmetic operators                          *
 ********************************************************/
inline Medium2D_State Medium2D_State::operator -(const Medium2D_State &U) const{
  Medium2D_State Temp(*this);
  Temp -= U;
  return Temp;
}

inline Medium2D_State Medium2D_State::operator +(const Medium2D_State &U) const{
  Medium2D_State Temp(*this);
  Temp += U;
  return Temp;
}

inline Medium2D_State Medium2D_State::operator *(const double &a) const{
  Medium2D_State Temp(*this);
  Temp *= a;
  return Temp;
}

inline Medium2D_State operator *(const double &a, const Medium2D_State &U) {
  Medium2D_State Temp(U);
  Temp *= a;
  return Temp;
}

inline Medium2D_State Medium2D_State::operator /(const double &a) const{
  Medium2D_State Temp(*this);
  Temp /= a;
  return Temp;
}


/********************************************************
 * Shortcut arithmetic operators                        *
 ********************************************************/
inline Medium2D_State& Medium2D_State::operator +=(const Medium2D_State &U){
  for(int i=0; i<Nband; i++) {
    kappa[i] += U.kappa[i];  
    sigma[i] += U.sigma[i];
    Ib[i]    += U.Ib[i];  
  }
  return (*this);
}

inline Medium2D_State& Medium2D_State::operator -=(const Medium2D_State &U) {
  for(int i=0; i<Nband; i++) {
    kappa[i] -= U.kappa[i];  
    sigma[i] -= U.sigma[i];
    Ib[i]    -= U.Ib[i];  
  }
  return (*this);
}

inline Medium2D_State& Medium2D_State::operator *=(const double &a) {
  for(int i=0; i<Nband; i++) {
    kappa[i] *= a;  
    sigma[i] *= a;
    Ib[i]    *= a;  
  }
  return (*this);
}

inline Medium2D_State& Medium2D_State::operator /=(const double &a) {
  for(int i=0; i<Nband; i++) {
    kappa[i] /= a;  
    sigma[i] /= a;
    Ib[i]    /= a;  
  }
  return (*this);
}

/********************************************************
 * Input/Output operators.                              *
 ********************************************************/
inline ostream& operator << (ostream &out_file, const Medium2D_State &U) 
{
  out_file.precision(10);
  out_file.setf(ios::scientific);
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.Ib[i];
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.kappa[i];
  for( int i=0; i<U.Nband; i++) out_file<<" "<<U.sigma[i];
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream& operator >> (istream &in_file,  Medium2D_State &U) 
{
  in_file.setf(ios::skipws);
  for( int i=0; i<U.Nband; i++) in_file >> U.Ib[i];
  for( int i=0; i<U.Nband; i++) in_file >> U.kappa[i];
  for( int i=0; i<U.Nband; i++) in_file >> U.sigma[i];
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif // _Medium2D_STATE_INCLUDED 
