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
 * Class Declaration                                    *
 ********************************************************/
class Medium2D_State;

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
#include "FieldData.h"

/********************************************************
 * Necessary Medium2D Specific Constants                *
 ********************************************************/

// Absorbsion model type 
enum Gas_Models { RTE2D_ABSORB_GRAY,
                  RTE2D_ABSORB_SNBCK };


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
  static FieldData* AbsorptionData;  //!< container class for absorbsion data
  static FieldData* ScatterData;     //!< container class for scatter data
  static FieldData* BlackbodyData;   //!< container class for blackbody data
  static TSpecificFunctor<FieldData>* func_kappa; //!< functor to return absorption coef
  static TSpecificFunctor<FieldData>* func_sigma; //!< functor to return scatter coef
  static TSpecificFunctor<FieldData>* func_Ib;    //!< functor to return blackbody
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
  //! memory allocation / deallocation for the field data arrays
  static void AllocateStatic();
  static void DeallocateStatic();
  //! memory allocation / deallocation for the SNBCK data object
  static void AllocateSNBCK();
  static void DeallocateSNBCK();
  //! deallocate all static variables
  static void DeallocateAllStatic() { DeallocateSNBCK(); DeallocateStatic(); }
  //@}

  //@{ @name State functions.
  //! Return extinction coefficient
  double beta(const int &v) const { return (kappa[v] + sigma[v]); }

  //! Return the state field
  double AbsorptionField(const Vector2D & r, const int &v)
  { return func_kappa[v](r); }
  double ScatterField(const Vector2D & r, const int &v)
  { return func_sigma[v](r); }
  double BlackbodyField(const Vector2D & r, const int &v)
  { return func_Ib[v](r); }

  //! Compute medium state at location
  friend void GetState(Medium2D_State &M, const Vector2D &r);
  void GetState(const Vector2D &r);
  //@}


  //@{ @name Static functions
  //! Setup function
  static void SetupStatic( const int &i_Absorb_Type, 
			   const SNBCK_Input_Parameters &SNBCK_IP,
			   const char* PATH);

  //! Set all all fields to the same function
  static void SetAllFieldsConstant(const Medium2D_State &M);
  //@}

  //@{ @name Index operator.
  double& operator[](int index);
  const double& operator[](int index) const;
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


inline void Medium2D_State :: AllocateStatic()
{
  // deallocate first
  DeallocateStatic();

  // allocate
  AbsorptionData = new FieldData[Nband];
  ScatterData    = new FieldData[Nband];
  BlackbodyData  = new FieldData[Nband];
  func_kappa     = new TSpecificFunctor<FieldData>[Nband];
  func_sigma     = new TSpecificFunctor<FieldData>[Nband];
  func_Ib        = new TSpecificFunctor<FieldData>[Nband];
}

inline void Medium2D_State :: DeallocateStatic()
{
  if ( AbsorptionData != NULL ) { delete[] AbsorptionData; AbsorptionData = NULL; }
  if ( ScatterData    != NULL ) { delete[] ScatterData;    ScatterData    = NULL; }
  if ( BlackbodyData  != NULL ) { delete[] BlackbodyData;  BlackbodyData  = NULL; }
  if ( func_kappa     != NULL ) { delete[] func_kappa;     func_kappa     = NULL; }
  if ( func_sigma     != NULL ) { delete[] func_sigma;     func_sigma     = NULL; }
  if ( func_Ib        != NULL ) { delete[] func_Ib;        func_Ib        = NULL; }
}

inline void Medium2D_State :: AllocateSNBCK()
{ 
  DeallocateSNBCK();
  SNBCKdata = new SNBCK; 
}

inline void Medium2D_State :: DeallocateSNBCK()
{ 
  if ( SNBCKdata != NULL ) { delete SNBCKdata; SNBCKdata = NULL;}  
}


/********************************************************
 * Compute values and initialize state.                 *
 ********************************************************/
inline void Medium2D_State :: SetInitialValues( const double &Pressure,
						const double &Temperature,
						const double &xco,
						const double &xh2o,
						const double &xco2,
						const double &xo2,
						const double &fsoot,
						const double &AbsorptionCoef,
						const double &ScatteringCoef)
{

  //------------------------------------------------
  // Absorbsion coefficient, Blackbody intensity 
  //------------------------------------------------
  // Use SNBCK
  if (Absorb_Type == RTE2D_ABSORB_SNBCK) {
    SNBCKdata->CalculateAbsorb( Pressure/PRESSURE_STDATM, //[atm]
				Temperature,              //[K]
				xco,
				xh2o,
				xco2,
				xo2,
				fsoot,
				kappa );
    SNBCKdata->CalculatePlanck( Temperature, Ib );

  // Use Gray Gas (ie. constant)
  } else if (Absorb_Type == RTE2D_ABSORB_GRAY) {
    for (int v=0; v<Nband; v++) {
      kappa[v] = AbsorptionCoef;
      Ib   [v] = BlackBody(Temperature);
    } // endfor

  // error
  } else{
    cerr << "Medium2D_State::SetInitialValues() - Invalid flag for Absorbsion model\n";
    exit(1);
  }

  //------------------------------------------------
  // Scattering coefficient 
  //------------------------------------------------
  // scattering coefficient always assumed gray
  for (int v=0; v<Nband; v++) { sigma[v] = ScatteringCoef; }

}



/********************************************************
 * Setup Static variables.                              *
 ********************************************************/
inline void Medium2D_State :: SetupStatic( const int &i_Absorb_Type, 
					   const SNBCK_Input_Parameters &SNBCK_IP,
					   const char* PATH) {

  //------------------------------------------------
  // Absorbsion model 
  //------------------------------------------------
  // set the absorption type flag
  Absorb_Type = i_Absorb_Type;
  
  // GRAY
  if (Absorb_Type == RTE2D_ABSORB_GRAY) {
    Nband = 1;

  // SNBCK
  } else if (Absorb_Type == RTE2D_ABSORB_SNBCK) {
    AllocateSNBCK();
    SNBCKdata->Setup(SNBCK_IP, PATH);
    Nband = SNBCKdata->NumVar();

  // ERROR
  } else {
    cerr << "Medium2D_State::SetupState - Invalid flag for gas type\n";
    exit(-1);
  } // endif

  // set number of variables
  NUM_VAR_MEDIUM2D = 3*Nband;

  //------------------------------------------------
  // Allocate Scalar fields 
  //------------------------------------------------
  Medium2D_State::AllocateStatic();

}


/********************************************************
 * Set all scalar fields to the same function.          *
 ********************************************************/
inline void Medium2D_State :: SetAllFieldsConstant(const Medium2D_State &M)
{
  // assign all the fields
  for (int i=0; i<Nband; i++) {

    // set constnat value
    AbsorptionData[i].val = M.kappa[i];
    ScatterData[i].val = M.sigma[i];
    BlackbodyData[i].val = M.Ib[i];

    // set field descriptor
    func_kappa[i].Set(&AbsorptionData[i], &FieldData::Constant);
    func_sigma[i].Set(&ScatterData[i], &FieldData::Constant);
    func_Ib[i].Set(&BlackbodyData[i], &FieldData::Constant);
  }
}



/********************************************************
 * Compute medium state at location.                    *
 ********************************************************/
inline void GetState(Medium2D_State &M, const Vector2D &r) 
{


  // compute
  for (int i=0; i<M.Nband; i++){
    M.kappa[i] = M.AbsorptionField(r,i);
    M.sigma[i] = M.ScatterField(r,i);
    M.Ib[i]    = M.BlackbodyField(r,i);
  } // endfor
}

inline void Medium2D_State :: GetState(const Vector2D &r) 
{
  for (int i=0; i<Nband; i++){
    kappa[i] = AbsorptionField(r,i);
    sigma[i] = ScatterField(r,i);
    Ib[i]    = BlackbodyField(r,i);
  } // endfor
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
