/********************* NASARP1311data.h **************************************
   This header file is for a class extracts species thermodynamic and transport
   information using the McBride & Gordon "NASA RP 1311" data.  The data is
   stored as polynomial functions of temperature with the coefficients in 
   the files "thermo.inp" and "trans.inp".

   This data corresponds to a "thermally perfect" gas which by definition
   (Anderson) is one where Cp and Cv are functions of temperature.

        data files 
               thermo.inp
               trans.inp

        assosiated class file: 
               NASARP1311class.cc
                
*****************************************************************************/

//SETUP
#ifndef _NASARP1311_DATA_INCLUDED
#define _NASARP1311_DATA_INCLUDED

#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>      //for stringstreams
#include <cstdlib> 	//for system command calls and getenv()
#include <cstdio>
#include <iomanip>

using namespace std;

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#define	INPUT_PARAMETER_LENGTH    128

/*************************************************************************
************** THEMOCOEF CLASS DEFINITION ********************************
  This is basically a data class that stores the thermodynamic data for 
  each temperatures range.
 
*************************************************************************
*************************************************************************/
class thermocoef{
  
private:
  
  //range of temperature for set of data
  double low_range;
  double high_range;
  
  //int numexponents;  //check, should be 7!
  double deltaHref;  //{H(298.15)-H(0)} J/mol
  
  // 0-6 thermo coefficients 
  // 7-8 b1, b2, integration constants
  double thermo_coef[9];

public:
  //default constructor
  //thermocoef::thermocoef()

  //data constructors input
  void Low_range_in(double L) {low_range = L;}
  void High_range_in(double H) {high_range = H;}
  void DeltaHref_in(double d)  {deltaHref = d;}
  void Thermo_coef_in(double *);
  
  //data constructors output
  double Low_range() {return low_range;}
  double High_range() {return high_range;}
  double DeltaHref() {return deltaHref;}
  //i = 0->8
  double Thermo_coef(int i);
  
  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const thermocoef &W);
  //  friend istream& operator >> (istream &in_file,  thermocoef &W);

  //destructor
  //thermocoef::~thermocoef();

};

/**************************************************************************
********************* THERMOCOEF CLASS CONSTRUCTORS ***********************
***************************************************************************/

/********** Read in Data *************************************************/
inline  void thermocoef::Thermo_coef_in(double *co){ 
  for( int i=0; i<9; i++){
    thermo_coef[i] = co[i];
    //cout<<co[i]<<endl;
  }
}

/*********** Ouput Data *************************************************/
inline double thermocoef::Thermo_coef(int i){
  //   if( i >= 0 && i <= 8){
  // shout check but too much overhead!
  return thermo_coef[i];
}

/*********** I/O Operator ***********************************************/
inline ostream &operator << (ostream &out_file, const thermocoef &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.low_range  << " " << W.high_range << " " 
	   << W.deltaHref<<endl;
  for( int i=0; i<9; i++){
    out_file<<" "<<W.thermo_coef[i];
  }
  out_file.unsetf(ios::scientific);
  return (out_file);
}

/*************************************************************************
************** TRANSCOEF CLASS DEFINITION ********************************
  This is basically a data class that stores the transport data for 
  each temperatures range.
 
**************************************************************************
*************************************************************************/
class transcoef{
  
 private:
  
  //range of temperature for set of data
  double low_range;
  double high_range;

  //int numexponents;  //check, should be 4!
 
  // 0-3 tansport coefficients 
  double trans_coef[4];

 public:
  //default constructor
  //transcoef::transcoef()

  //data constructors input
  void Low_range_in(double L) {low_range = L;}
  void High_range_in(double H) {high_range = H;}
  void Trans_coef_in(double *); 
  
  //data constructors output
  double Low_range()  {return low_range;}
  double High_range() {return high_range;}
  
  //int 1->4
  double Trans_coef(int i);
   
  /* Input-output operators. */ 
  friend ostream& operator << (ostream &out_file, const transcoef &W);
  
  //destructor
  //transcoef::~transcoef();

};

/**************************************************************************
********************* TRANSCOEF CLASS CONSTRUCTORS ***********************
***************************************************************************/

inline ostream &operator << (ostream &out_file, const transcoef &W) {
  out_file.setf(ios::scientific);
  out_file << " " << W.low_range  << " " << W.high_range <<"   "; 
    //	   << W.numexponents <<endl;
  for( int i=0; i<4; i++){
    out_file<<" "<<W.trans_coef[i];
  }
  out_file.unsetf(ios::scientific);
  return (out_file);
}

/********** Read in Data **************************************************/
inline  void transcoef::Trans_coef_in(double *co){ 
  for( int i=0; i<4; i++){
    trans_coef[i] = co[i];
  }
}

/*********** Ouput Data *************************************************/
inline double transcoef::Trans_coef(int i){
  //   if( i >= 0 && i <= 8){
  // shout check but too much overhead! 
  return trans_coef[i];
}


/*************************************************************************
************** NASAPR1311data CLASS DEFINITION ***************************
  This class parses the species thermodynamic and transport data files
  and provides constructors to retrieve the following physical data.
  
  Heat of formation:   deltaHf   
  Molecular Mass:          M    kg/mol
  Species Gas Constant:    R    J/(kg*K)
  {H(298.15)-H(0)}              J/kg
  
  Temperature(Kelvin) Dependent data: f(T) 
   
  Heat Capacity (const Pressure):  Cp  J/(kg*K)
  Heat Capacity (const Volume):    Cv  J/(kg*K)
  Specific Heat Ratio:             gamma 
  Specific Enthalpy:               H   J/kg
  Specific Entropy:                S   J/(kg*K)
  Specific Internal Energy:        U   J/kg
  Viscosity:                       mu  kg/(m*s) N*s/m^2  
  Thermal Conductivity:            k   N/(s*K)  W.(m*K)

  It uses the thermocoef & transcoef classes to handle multiple data sets 
  for the same species over different temperature ranges.

  If the temperature range is below the lowest data set bounds
  the low_range values are used as an approximation.

  Be careful with the h and s values, its based on Tref=298.15K so in
  order to match absolute reference data you need to add the DeltaHref
  amount.  Also the direct specific heats are "Total" in the sense 
  they contain the heat of formation in them already.  So to 
  calculate the "sensible" specific enthalpies and "sensible" specific
  internal energy use the following relations.

  h = H + DeltaHref - heatofform
  u = U + DeltaHref - heatofform

  Plot_data uses gnuplot to plot the physical data polynomials for
  all the temperature ranges available.  The "types" available
  are Cp,H,S,mu, and k. The output is a postscript file.

************************************************************************* 
*************************************************************************/ 
class NASARP1311data{
 
  private: 

  //data file name and path 
  //const char datafilename[15] = "thermo.inp";
  char datafilename_thermo[INPUT_PARAMETER_LENGTH]; 
  char datafilename_trans[INPUT_PARAMETER_LENGTH];

  // species name
  string species;
  double mol_mass;
  double heatofform;
  
  //number of temperature ranges
  int temp_ranges_thermo; 
  int temp_ranges_V;
  int temp_ranges_C;
 
  //sets of thermodynamic coefficients 
  //for each temperature range
  thermocoef *thermo_data;
  
  //sets of transport properties for
  //each temperature range
  transcoef *trans_thermconduct;
  transcoef *trans_viscosity;

  //static NASARP1311data default_data;
  public:
  //default constructor
  NASARP1311data(){ 
    strcpy(datafilename_thermo,"thermo.inp"); 
    strcpy(datafilename_trans,"trans.inp");
    species = "O2";
    mol_mass = 32.0;
    heatofform = 0.0;
    temp_ranges_thermo=1;
    temp_ranges_V=1;
    temp_ranges_C=1; 
    thermo_data=NULL; 
    trans_thermconduct=NULL; 
    trans_viscosity=NULL;   
  } 
  
  //read species parameters from data file
  void Getdata(string, const char *);
  void Set_Path_Names(const char *);

  //Retrieve specfic data
  string Speciesname() {return species;}
  double Mol_mass() {return mol_mass;}
  double Heatofform() {return heatofform;}
  double Rs(); 
  int Temp_ranges_thermo() { return temp_ranges_thermo;}
  int Temp_ranges_V() { return temp_ranges_V;}
  int Temp_ranges_C() { return temp_ranges_C;}
  double DeltaHref() { return thermo_data[0].DeltaHref();}

  //find which set of coefficients to use (internal use)
  int Which_coef(double& Temp);
  int Which_coef(double& Temp, transcoef *trans, int temp_ranges);
  double Low_range(){ return thermo_data[0].Low_range(); }
  double High_range(){ return thermo_data[temp_ranges_thermo -1].High_range();}

  //find thermo data useing polynomials
  double HeatCapacity_p(double Temp);  //Cp 
  double Enthalpy(double Temp);        //H (J/mol)
  double Enthalpy_mol(double Temp);    //H (J/mol)
  double Enthalpy_prime(double Temp);  // dh/dT
  double Entropy(double Temp);         //S (J/kgK) 
  double Entropy_mol(double Temp);     //S (J/molK)
  
  //find transport propteries using coefficients
  double Viscosity(double Temp);
  double ThermalConduct(double Temp);
  double dViscositydT(double Temp);
  double dKappadT(double Temp);

  //furthur constructors based on data
  double InternalEnergy(double Temp);  //U
  double HeatCapacity_v(double Temp);  //Cp 
  double HeatRatio(double Temp);       //gamma 
  double Prandtl(double Temp);         //Pr

  // Input-output operators. 
  friend ostream& operator << (ostream &out_file, const NASARP1311data &W);
  //  friend istream& operator >> (istream &in_file,  thermocoef &W);
  
  //plot out curve fits
  //void Plot_data(string type);

  void deallocate();
  //Destructor 
  ~NASARP1311data(){ deallocate();}
};


/**************************************************************************
********************* NASARP1311data CLASS CONSTRUCTORS *******************
***************************************************************************/

/**************** Data file Path Names ***************************/
inline void NASARP1311data::Set_Path_Names(const char *CFFC_path){
  
  //Set NASA data path 
  strcpy(datafilename_thermo, CFFC_path);
  strcat(datafilename_thermo, "/data/NASARP1311/thermo.inp");
  strcpy(datafilename_trans, CFFC_path);
  strcat(datafilename_trans, "/data/NASARP1311/trans.inp");

}

/********************* Species Gas Constant ***********************/
inline double NASARP1311data::Rs(){ 
  return R_UNIVERSAL/mol_mass;
}

/********************** Internal energy ***************************/
inline double NASARP1311data::InternalEnergy(double Temp){
  //internal energy
  // u = h - P/rho
  return (Enthalpy(Temp) - (R_UNIVERSAL/mol_mass)*Temp);  
}

inline double NASARP1311data::HeatCapacity_v(double Temp){ 
  //Cv = Cp -R;
  return HeatCapacity_p(Temp) - Rs();
} 

inline double NASARP1311data::HeatRatio(double Temp){
  //gamma = Cv/Cp
  return HeatCapacity_p(Temp)/HeatCapacity_v(Temp);
}

inline double NASARP1311data::Prandtl(double Temp){
  //Pr = Cp*mu/k
  return HeatCapacity_p(Temp)*Viscosity(Temp)/ThermalConduct(Temp);
}

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
/**************************** DESTRUCTOR ****************************************/
inline void NASARP1311data::deallocate()
{ 
  //from Getdata
  if( thermo_data !=NULL){ 
    delete[] thermo_data;
  }
  if( trans_viscosity !=NULL){
    delete[] trans_viscosity; 
  }
  if( trans_thermconduct !=NULL){
    delete[] trans_thermconduct;
  }
}

/********************** I/O Operator *****************************/
inline ostream &operator << (ostream &out_file, const NASARP1311data &W) {
  out_file.setf(ios::scientific);
  out_file << "\n " << W.species << " " << W.heatofform << " " 
	   << W.mol_mass << " " << W.temp_ranges_thermo<<endl;
  
  //Thermodynamic
  cout<<"\n Thermodynamic";
  for( int i=0; i<W.temp_ranges_thermo; i++){
    out_file<<endl<<W.thermo_data[i];
  }
   
  //using string to double from standard C library
  //const char* Temp = s.c_str();  //Transport -> viscosity 
  cout<<"\n Viscosity";
  for( int i=0; i<W.temp_ranges_V; i++){
    out_file<<endl<<W.trans_viscosity[i];
  }

  //Transport -> thermal conductivity
  cout<<"\n Thermal Conducivity ";
  for( int i=0; i<W.temp_ranges_C; i++){
    out_file<<endl<<W.trans_thermconduct[i];
  }
  
  out_file.unsetf(ios::scientific);
  return (out_file);
}//end of dataclass constructors

#endif /* END of _NASARP1311_DATA_INCLUDED */

//conversion functions
extern int stoi(const string);
extern double stof(const string);
