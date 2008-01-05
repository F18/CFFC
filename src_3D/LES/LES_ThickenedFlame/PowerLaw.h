/****************** PowerLaw.h **************************************
  This class defines the subfilter scale model for the turbulent flame
  speed. The model is based on writing the unresolved flame surface
  density in terms of a general power-law expression that involves an
  inner and an outer cutoff scales.  
********************************************************************/

#ifndef _POWER_LAW_INCLUDED
#define _POWER_LAW_INCLUDED

#include <iostream>
#include <fstream>

using namespace std;

// Required CFFC header files
#ifndef _MATH_MACROS_INCLUDED
#include "../../Math/Math.h"
#endif //MATH_MACROS_INCLUDED

#ifndef _GAS_CONSTANTS
#include "../../Physics/GasConstants.h"
#endif


//! Model constant for the SFS root mean square velocity. 
const double C_sfs_rms = 2.0;      

//! Factors used in the dynamic thickening formulation.
const double Beta_factor = 2.0;
const double Gamma_factor = 0.5;




/************************************************************************
********************* POWER_LAW CLASS DEFINITION ************************

  Wrinkling factor:                 WF          
  SFS turbulent flame speed:        sfs_s_turb        m/s
  Filter size:                      filter_size        m
  SFS Reynolds number:              sfs_Re      
  SFS turbulent speed:              u_prime           m/s
  Mean curvature:                   mean_curvature    1/m
  Viscosity:                        mu              kg/(m*s)
  Density:                          rho              kg/m^3
  Laminar flame speed:              lam_speed         m/s
  Laminar flame thickness:          lam_thickness     m/s

**************************************************************************
**************************************************************************/

class PowerLaw{
private:
public:
  
  double WF;             // SFS wrinkling factor
  double TF;             // Used for dynamic thickening

  static const int n = 8;  // Number of cells across the flame front                
  
  // constructors
  PowerLaw() : WF(ONE), TF(ONE) { }
  PowerLaw(const double &wf, const double &tf) : WF(wf), TF(tf) { }


  void unphysical_check(const double &TFactor);

  /********************  Filter size  *********************************  
     This is the size of the filter used for the thickened flame,
     which is greater than, and different from, the LES filter.
  ********************************************************************/  
  double filter_size(const double &lam_thickness) const;
  double mesh_size(const double &cell_volume) const;  


  /*******************  Thickening factor ****************************/
  void thickening_factor(const double &cell_volume,
			 const double &lam_thickness);

  void TF_function(const double &TFactor,
		   const double &progress_variable);

  void TF_function(const double &fuel_mass_fraction,
		   const double &temperature, 
		   const double &TFactor);

  void TF_function(const double &fuel_mass_fraction,
		   const double &oxi_mass_fraction,
		   const double &temperature, 
		   const double &TFactor);

  double sensor_omega(const double &fuel_mass_fraction,
		      const double &temperature);

  double sensor_omega(const double &fuel_mass_fraction,
		      const double &oxi_mass_fraction,
		      const double &temperature);


  /***************** SFS turbulent velocity **************************/
  double u_prime(const double &cell_size,  //const double &lam_thickness, 
                 const double &lapl_vor) const;
 

  /***************** Efficiency function *****************************
     Efficiency function for the effective straining that accounts 
     for the effects of all relevant scales.
  *******************************************************************/ 
  double efficiency_function(const double &lam_speed,
                             const double &lam_thickness,
                             const double &cell_size,
                             const double &lapl_vor,
                             const double &rho,
                             const double &mu) const;
  

  /*****************    Mean curvature   ****************************
    This function determines the mean curvature, defined as the 
    inverse of the inner cutoff scale.
  *******************************************************************/
  double mean_curvature(const double &lam_speed,
			const double &lam_thickness,
			const double &cell_size,
			const double &lapl_vor,
			const double &rho,
			const double &mu) const;


  /****************  Wrinkling factor  ******************************
    Function that calculates the wrinkling factor, which is the 
    ratio of the turbulent flame speed to the laminar one. It is used 
    as the efficiency factor in the thickened flame model.
  *******************************************************************/  
  void wrinkling_factor(const double &lam_speed,
			const double &lam_thickness,
			const double &cell_size,
			const double &lapl_vor,
			const double &rho,
			const double &mu);  


  /***************  Turbulent flame speed   *************************
     Function that calculates the SFS turbulent flame speed.
  *******************************************************************/
  double turb_speed(const double &lam_speed);
  
  
  
  /****************** Operator Overloading **************************/
  PowerLaw operator +(const PowerLaw &PL) const;
  PowerLaw operator -(const PowerLaw &PL) const;
  double  operator *(const PowerLaw &PL) const;
  PowerLaw operator *(const double &a) const;
  friend  PowerLaw operator *(const double &a, const PowerLaw &PL);
  PowerLaw operator /( const double &a) const;
  

  /* Unary arithmetic operators. */
  PowerLaw operator +(void) const;
  PowerLaw operator -(void) const;
  
  /* Assignment Operator */
  PowerLaw &operator =(const PowerLaw &PL); 

  /* Shortcut arithmetic operators. */
  PowerLaw &operator +=(const PowerLaw &PL);
  PowerLaw &operator -=(const PowerLaw &PL);
  
  /* Relational operators. */
  int operator ==(const PowerLaw &PL) const;
  int operator !=(const PowerLaw &PL) const; 
  
  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const PowerLaw &PL);
  friend istream& operator >> (istream &in_file,  PowerLaw &PL);

  // Use default destructor
   
};




/*----------------------------------------------------------*\
          Inline member functions
\*----------------------------------------------------------*/

inline void PowerLaw::unphysical_check(const double &TFactor) {
  if (TF < ONE-NANO) { cout << "\n WARNING: TF smaller than allowed = " << TF << flush; TF = ONE; }
  if (TF > TFactor+NANO) { cout << "\n WARNING: TF greater than allowed = " << TF << flush; TF = TFactor; }
  if (WF < ONE-NANO) { cout << "\n WARNING: WF smaller than allowed = " << WF << flush; WF = ONE; }
  if (WF > pow(TFactor, 2.0/3.0)) { cout << "\n WARNING: WF greater than allowed = " << WF << flush; WF = pow(TFactor, 2.0/3.0); }
}


//Thickening factor in terms of the mesh size and number of cells across the flame front
inline  void PowerLaw::thickening_factor(const double &cell_volume,
					 const double &lam_thickness) {
  
  TF = ceil(double(n)*mesh_size(cell_volume)/lam_thickness);
  if (TF >= TEN) cout << "\n WARNING: TF = " << TF << " greater than 10." << flush;
  //cout << "\n In PowerLaw n = " << n << "\t TF = " << TF << flush;
}



inline void  PowerLaw::TF_function(const double &TFactor,
				   const double &progress_variable) {

    double z = (progress_variable - 0.5)/(0.14*sqrt(TWO));
    TF = (TFactor-ONE)*exp(-z*z) + ONE ;
}

inline void  PowerLaw::TF_function(const double &fuel_mass_fraction,
				   const double &temperature, 
				   const double &TFactor) {

  double sensor_max, sensor_value;
  
  // Maximum value of the sensor function
  // CH4, equivalence ratio phi = 0.6    =>  932.615
  // CH4, equivalence ratio phi = 0.8    =>  
  // CH4, equivalence ratio phi = 1.0    =>  

  sensor_max = 932.615;   //185935.0;    //7.31337; //  6.22006;  
  sensor_value = sensor_omega(fuel_mass_fraction, temperature);

  TF = (TFactor-ONE)*erf(Beta_factor*sensor_value/sensor_max) + ONE;
}


// overloaded for fuel and oxidizer
inline void  PowerLaw::TF_function(const double &fuel_mass_fraction,
				   const double &oxi_mass_fraction,
				   const double &temperature, 
				   const double &TFactor) {

  double sensor_max, sensor_value;

  // Maximum value of the sensor function
  // CH4, equivalence ratio phi = 0.6    =>  103.304 
  // CH4, equivalence ratio phi = 0.8    =>  394.246
  // CH4, equivalence ratio phi = 1.0    =>  886.517
 
  sensor_max = 8.86449e-07;    //25174.0;    //7.31337; //  6.22006;  
  sensor_value = sensor_omega(fuel_mass_fraction, oxi_mass_fraction, temperature);

  TF = (TFactor-ONE)*erf(Beta_factor*sensor_value/sensor_max) + ONE;
}


inline double PowerLaw::sensor_omega(const double &fuel_mass_fraction,
				     const double &temperature) {
  

  return (fuel_mass_fraction/*NANO*/)*exp(-Gamma_factor*48400.0*4.1868/(R_UNIVERSAL*temperature));
}


// overloaded for fuel and oxidizer
inline double PowerLaw::sensor_omega(const double &fuel_mass_fraction,
				     const double &oxi_mass_fraction,
				     const double &temperature) {

  return (fuel_mass_fraction*oxi_mass_fraction/*NANO*/)*exp(-Gamma_factor*48400.0*4.1868/(R_UNIVERSAL*temperature));
}


// Mesh size
inline double PowerLaw::mesh_size(const double &cell_volume) const {
  return pow(cell_volume, 1.0/3.0);
}


// Filter size
inline  double PowerLaw::filter_size(const double &lam_thickness) const {

  //const int np = 10;
  //return  (double(np)*filter_width);
  if (TF < ONE) { cout << "\nWarning: Thickening factor less than 1." << flush; }
  return (TF*lam_thickness);
}


// SFS turbulence intensity
inline double PowerLaw::u_prime(const double &cell_size,  //const double &lam_thickness,    
				const double &lapl_vor) const {
  double fs; 
  fs = cell_size;  //fs = filter_size(lam_thickness);        
  

  return  (C_sfs_rms*fs*fs*fs*lapl_vor);  //return  0.5;
}


// SFS turbulent flame speed
inline double PowerLaw::turb_speed(const double &lam_speed) {
  return (WF*lam_speed);
}




/*----------------------------------------------------------*\
                        Operators
\*----------------------------------------------------------*/

//----------------- Addition -----------------------------//
inline PowerLaw PowerLaw::operator +(const PowerLaw &PL) const{
  return(PowerLaw(PL.WF+WF, PL.TF+TF));
}

//------------------ Subtraction -------------------------//
inline PowerLaw PowerLaw::operator -(const PowerLaw &PL) const{
    return(PowerLaw(WF-PL.WF, TF-PL.TF));
}

//----------------- Inner Product ------------------------//
inline double PowerLaw::operator *(const PowerLaw &PL) const{
      return(WF*PL.WF + TF*PL.TF); 
}

//---------------- Scalar Multiplication -----------------//
inline PowerLaw PowerLaw::operator *(const double &a) const{
  return(PowerLaw(a*WF, a*TF));
}

inline PowerLaw operator *(const double &a, const PowerLaw &PL){
  return(PowerLaw(a*PL.WF, a*PL.TF));
}
 
//--------------- Scalar Division ------------------------//
inline PowerLaw PowerLaw::operator /(const double &a) const{
  return(PowerLaw(WF/a, TF/a));
}
 
//------------- Unary arithmetic operators ----------------//
inline PowerLaw PowerLaw::operator +(void) const{
  return(PowerLaw(WF, TF)); 
}
 
inline PowerLaw PowerLaw::operator -(void) const{
  return(PowerLaw( -WF, -TF)); 
}
 
//----------------- Assignment ----------------------------//
inline PowerLaw &PowerLaw::operator =(const PowerLaw &PL){ 
  WF = PL.WF;
  TF = PL.TF; 
  return(*this);
}

//----------- Shortcut arithmetic operators ---------------//
inline PowerLaw &PowerLaw::operator +=(const PowerLaw &PL){
  WF += PL.WF;
  TF += PL.TF;
  return(*this);
}

inline PowerLaw &PowerLaw::operator -=(const PowerLaw &PL){
  WF -= PL.WF;
  TF -= PL.TF;
  return(*this);
}

//----------------- Relational operators ------------------// 
inline int PowerLaw::operator ==(const PowerLaw &PL) const{
  return (WF == PL.WF && TF == PL.TF);
}

inline int PowerLaw::operator !=(const PowerLaw &PL) const{
   return (WF != PL.WF || TF != PL.TF);
}


//------------------ I/O operators ------------------------//
inline ostream &operator << (ostream &out_file, const PowerLaw &PL) {
  out_file.setf(ios::scientific);
  out_file <<" "<< PL.WF <<" "<< PL.TF;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, PowerLaw &PL) {
  in_file.setf(ios::skipws);
  in_file >> PL.WF >> PL.TF;
  in_file.unsetf(ios::skipws);
  return (in_file);
}


#endif // POWER_LAW_INCLUDED




