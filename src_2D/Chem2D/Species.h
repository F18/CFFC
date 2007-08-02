/************************* Species.h **************************************
  This header file is for the Species class as part of Chem2D. Essentially
  it is used by Chem2DState classes to hold the species specific data, which
  so far is only the mass fraction, c.  The global thermodynamic and 
  transport data is stored statically for each species using the NASARP1311data 
  class.
               
*****************************************************************************/

#ifndef _SPECIES_INCLUDED   
#define _SPECIES_INCLUDED

#include <cmath>
#include <iostream>

using namespace std;

#ifndef _VECTOR2D_INCLUDED
#include "../Math/Vector2D.h"
#endif //_VECTOR2D_INCLUDED

/*************************************************************************
 ******************** SPECIES CLASS DEFINITION ***************************
  DATA:
   -> c                   mass fraction of species
   -> diffusion_coef      diffusion coefficient
   -> gradc               gradients of mass fraction (Vector)

  CONSTRUCTORS, etc.. :
   -> operator overloads
************************************************************************* 
*************************************************************************/ 
class Species {
 
 private: 
 public: 
  //mass fraction
  double c;
  double diffusion_coef;
  Vector2D gradc;
 
  //read species parameters from data file
  Species(){c=ONE; diffusion_coef=ZERO; gradc.x=ZERO; gradc.y=ZERO;}
  Species(const double &frac){c=frac; diffusion_coef=ZERO; gradc.x=ZERO; gradc.y=ZERO;}
  Species(const double &frac, const double &dcoef, const Vector2D &gc)
  {c=frac; diffusion_coef=dcoef; gradc=gc; }
  
  /*************** VACUUM OPERATOR *********************/
  void Vacuum(){c=ZERO; diffusion_coef=ZERO; gradc.zero();}

  /****************** Operator Overloading **************************/
  Species operator +(const Species &A) const;
  Species operator -(const Species &A) const;
  double  operator *(const Species &A) const;
  Species operator *(const double &a) const;
  friend  Species operator *(const double &a, const Species &A);
  Species operator /( const double &a) const;
  

  /* Unary arithmetic operators. */
  Species operator +(void) const;
  Species operator -(void) const;
  
  /* Assignment Operator */
  Species &operator =(const Species &A); 

  /* Shortcut arithmetic operators. */
  Species &operator +=(Species &A);
  Species &operator -=(Species &A);
  Species &operator *=(const double &a);
  Species &operator /=(const double &a);

  /* Relational operators. */
  int operator ==(const Species &A) const;
  int operator !=(const Species &A) const;
  
  /* Input-output operators. */
  friend ostream &operator << (ostream &out_file, const Species &B);
  friend istream &operator >> (istream &in_file,  Species &B);
  //Destructor 
  //Species::~Species();	   
};

/**************************************************************************
********************* Species CLASS CONSTRUCTORS **************************
**************************************************************************/

//----------------- Addition -----------------------------//
inline Species Species::operator +(const Species &A) const{
  return(Species(A.c+c));
// 		 A.diffusion_coef+diffusion_coef,
// 		 A.gradc+gradc));
}

//------------------ Subtraction -------------------------//
inline Species Species::operator -(const Species &A) const{
    return(Species(c-A.c));
// 		 diffusion_coef-A.diffusion_coef,
// 		 gradc-A.gradc));
}

//----------------- Inner Product ------------------------//
inline double Species::operator *(const Species &A) const{
      return(c*A.c); 
// 	 diffusion_coef*A.diffusion_coef +
// 	 gradc*A.gradc); 
}

//---------------- Scalar Multiplication -----------------//
inline Species Species::operator *(const double &a) const{
  return(Species(a*c));
// 		 a*diffusion_coef,
// 		 a*gradc));
}

inline Species operator *(const double &a, const Species &A){
  return(Species(a*A.c));
  // 		 a*A.diffusion_coef,
  // 		 a*A.gradc));
}
 
//--------------- Scalar Division ------------------------//
inline Species Species::operator /(const double &a) const{
  return(Species(c/a));
  // 		 diffusion_coef/a,
  // 		 gradc/a));
}
 
//------------- Unary arithmetic operators ----------------//
inline Species Species::operator +(void) const{
  return(Species(c)); //,diffusion_coef,gradc));
}
 
inline Species Species::operator -(void) const{
  return(Species(-c)); //,-diffusion_coef,-gradc));
}
 
//----------------- Assignment ----------------------------//
inline Species &Species::operator =(const Species &A){
  c = A.c; 
  diffusion_coef= A.diffusion_coef;
  gradc = A.gradc; 
  return(*this);
}

//----------- Shortcut arithmetic operators ---------------//
inline Species &Species::operator +=(Species &A){
  c += A.c; 
//   diffusion_coef += A.diffusion_coef;
//   gradc += A.gradc; 
  return(*this);
}

inline Species &Species::operator -=(Species &A){
  c -= A.c; 
//   diffusion_coef -= A.diffusion_coef;
//   gradc -= A.gradc; 
  return(*this);
}

inline Species& Species::operator *=(const double &a) {
  c *= a; 
  return *this;
}

inline Species& Species::operator /=(const double &a) {
  c /= a; 
  return *this;
}

//----------------- Relational operators ------------------// 
inline int Species::operator ==(const Species &A) const{
  return (c == A.c && diffusion_coef == A.diffusion_coef && gradc == A.gradc);
}

inline int Species::operator !=(const Species &A) const{
   return (c != A.c || diffusion_coef != A.diffusion_coef || gradc != A.gradc);
}


//------------------ Input-output operators ---------------// 
inline ostream &operator << (ostream &out_file, const Species &B){
  out_file.setf(ios::scientific);
  out_file<<" "<<B.c; //<<" "<<B.gradc<<" "<<B.diffusion_coef;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file,  Species &B){
  in_file.setf(ios::skipws);
  in_file>>B.c; //>>B.gradc>>B.diffusion_coef;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

#endif /* END of _SPECIES_INCLUDED */



