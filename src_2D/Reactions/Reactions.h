/****************** Reactions.h **************************************
  This class defines the Reaction mechanisms for the 
  2D Axisymmertric Navier-Stokes with Multiple Species solution.

  TODO:   - kb constructor's
          - user defind reaction set constructor
          - use int's instead of voids for error_flagging

***********************************************************************/
#ifndef _REACTIONS_INCLUDED
#define _REACTIONS_INCLUDED 

class React_data;
class Reaction_set;

// Required C++ libraries
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

#include "../Math/Math.h"
#include "../Math/Matrix.h"
#include "../Physics/GasConstants.h"
#include "../CFD/CFD.h"
#include "../Math/Complexify.h"
#include "15step.h"

// Cantera libraries
#ifdef _CANTERA_VERSION
#include <cantera/Cantera.h>      // main include
#include <cantera/IdealGasMix.h>  // reacting, ideal gas mixture class
#endif //_CANTERA_VERSION

//Calorie to Joule conversion
#define CAL_TO_JOULE 4.1868 // 4.1868 Joules = 1 calorie

/************************************************************************
 ***************** REACTION DATA CLASS **********************************
  This class is another data container used for the Arhenius coefficients
  used in caluclated the forward (kf) and backward (kb) reaction rates.

  kf = A(T^n)exp( -E/RT)

  kov = A(T^n)exp( -E/RT)[Fuel]^a[Oxidizer]^b

  Be careful of units, need to get k in (m^3/(mol*s)) for 2 body.
  and (m^6/(mol*s)) for 3rd body.  Most A coefficients are in terms
  of cm^3/(mol*s) so be sure to check to units.  Also the 
  if species concentrations are used, as in the kov, they are 
  are dimensionalized and that is why in the kov member functions
  below they are divided by 1e6.  Confusing yes, but that is the
  way that data is presnted.
 
    A : m^3/(mol*s) FYI 1e6cm^3 = 1m^3
    Ea: J/mol

************************************************************************
************************************************************************/
class React_data {
private: 
  string react; //reaction name

  double A;  //three Arhenius coefficients
  double Ab;
  double n;
  double E;
  
  //for keq formulation
  double nu_coef; //stoichmetric coefficients for Keq => Kc

protected:
public:
  //default constructor
  // name, A ,b ,E
  void set_data(string nam,double x,double y ,double z, int nu){
    react=nam, A=x; n=y; E=z; nu_coef=nu;}
  
  void set_data(string nam,double x, double xx, double y, double z, int nu)
  { react=nam; A=x; Ab=xx; n=y; E=z; nu_coef=nu;}
  
  //get reaction rate coefficients
  double kf(const double &Temp) const;
  double dkf_dT(const double &Temp) const;
  double kf(const double &Temp, const double &H2, const  double &O2, const double &N2)const; //for H2&O2
  double kb(const double &Temp) const;
  double dkb_dT(const double &Temp) const;

  template<class SOLN_pSTATE>  
  double keq(const SOLN_pSTATE &W, const double &Temp) const;

  //Determine change in Gibbs free energy
  template<class SOLN_pSTATE>  
  double deltaG(const SOLN_pSTATE &W, const double& Temp) const;

  string react_name()const{ return react;}

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const React_data &W);
  friend istream& operator >> (istream &in_file,  React_data &W);
};

/*************** forward reaction coef (kf) ****************************/
inline double React_data::kf(const double &Temp) const{
  return A*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/****** derivative of forward reaction coef (kf) wrt to Temperature *********/
inline double React_data::dkf_dT(const double &Temp) const{
  return A*pow(Temp,n-1)*exp(-E/(R_UNIVERSAL*Temp))*(n + E/(R_UNIVERSAL*Temp));
}

/*************** forward reaction kf for 2Step H2 ***********************/
inline double React_data::kf(const double &Temp, const double &H2, const double &O2, const double &N2) const{
  //equivalence ratio
  double stoich = 4.0/(32.0+3.76*28.0);
  double phi = (H2/(O2+N2))/stoich;
  if( phi < TOLER){
    phi = TOLER;
  }
  double Astar=ZERO;

  //using A for units conversion cm^3/mol*s -> m^3/mol*s;
  if(react == "H2O2_2step_1"){                          
    Astar = A*(8.917*phi + (31.433/phi) - 28.950)*1e47; 
  } else if (react == "H2O2_2step_2"){ 
    Astar = A*(2.0 + (1.333/phi) - 0.833*phi)*1e64; 
  }

//   cout<<endl<<phi<<" "<<A<<" "<<(-E/(R_UNIVERSAL*Temp))
//       <<" "<<R_UNIVERSAL<<" "<<E<<" "<<Temp<<" "<<Astar*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));

  return Astar*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/*************** backward reaction coef (kb) ****************************/
inline double React_data::kb(const double &Temp) const{
  return Ab*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/****** derivative of backward reaction coef (kb) wrt to Temperature *********/
inline double React_data::dkb_dT(const double &Temp) const{
   return Ab*pow(Temp,n-1)*exp(-E/(R_UNIVERSAL*Temp))*(n + E/(R_UNIVERSAL*Temp));
}

/***************** equilibrium coef Keq ********************************/
template<class SOLN_pSTATE>  
inline double React_data::keq(const SOLN_pSTATE &W, const double& Temp) const{
  // nu_coef is the stoichiometric coef. sum (Eqn 13.25, 13.26 Anderson)
  // Kp or Keq has units of Pressure Pa( N/m^2) so need to change to 
  // cgs units i.e *1e6
  return pow(R_UNIVERSAL*Temp, nu_coef)*exp(-deltaG(W, Temp)/(R_UNIVERSAL*Temp))*1e6;
}

/**********************************************************************
   Determine the Gibbs free energy change 
   deltaG = Sum(products) - Sum(reactants)
   
   Units:  J/mol
***********************************************************************/
template<class SOLN_pSTATE>  
inline double React_data::deltaG(const SOLN_pSTATE &W, const double& Temp) const{

  //need a general and more efficient system here !!!! SN
  // deltaG = sum (( product Gs) - (react Gs))

  //2STEP
  if(react == "H2O2_2step_1"){
    return TWO*W.specdata[2].GibbsFree(Temp) - ( W.specdata[0].GibbsFree(Temp) + W.specdata[1].GibbsFree(Temp));
  } else if (react == "H2O2_2step_2"){
    return TWO*W.specdata[3].GibbsFree(Temp) - ( TWO*W.specdata[2].GibbsFree(Temp) + W.specdata[0].GibbsFree(Temp));
 
    //8STEP
  } else if (react == "H2O2_8step_1"){
    return W.specdata[5].GibbsFree(Temp) + W.specdata[0].GibbsFree(Temp) - ( W.specdata[2].GibbsFree(Temp) + W.specdata[1].GibbsFree(Temp));
  } else if (react == "H2O2_8step_2"){
    return W.specdata[5].GibbsFree(Temp) + W.specdata[2].GibbsFree(Temp) - ( W.specdata[0].GibbsFree(Temp) + W.specdata[3].GibbsFree(Temp));
  } else if (react == "H2O2_8step_3"){
    return W.specdata[4].GibbsFree(Temp) + W.specdata[2].GibbsFree(Temp) - ( W.specdata[5].GibbsFree(Temp) + W.specdata[3].GibbsFree(Temp));
  } else if (react == "H2O2_8step_4"){
    return W.specdata[4].GibbsFree(Temp) + W.specdata[0].GibbsFree(Temp) - ( W.specdata[5].GibbsFree(Temp) + W.specdata[5].GibbsFree(Temp));
  } else if (react == "H2O2_8step_5"){
    return W.specdata[3].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp) - ( TWO*W.specdata[2].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp));
  } else if (react == "H2O2_8step_6"){
    return W.specdata[4].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp) - ( W.specdata[2].GibbsFree(Temp) + W.specdata[5].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp));
  } else if (react == "H2O2_8step_7"){
    return W.specdata[5].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp) - ( W.specdata[2].GibbsFree(Temp) + W.specdata[0].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp));
  } else if (react == "H2O2_8step_8"){
    return W.specdata[1].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp) - ( TWO*W.specdata[0].GibbsFree(Temp) + W.specdata[4].GibbsFree(Temp));
  } else {
    cerr<<" \n Missing deltaGibbs for "<<react;
    exit(1); 
  }
	     
}

/**************** I/O Operators ***************************************/
inline ostream &operator << (ostream &out_file, const React_data &W) {
  out_file.setf(ios::scientific);
  out_file <<"\n "<<W.react<<" A "<<W.A<<" n "<<W.n<<" E "<<W.E;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file, React_data &W) {
  in_file.setf(ios::skipws);
  in_file >> W.react >> W.A >> W.n >> W.E;
  in_file.unsetf(ios::skipws);
  return (in_file);
}


/**************************************************************************
********************* REACTIONS CLASS DEFINTION ***************************
 This class uses the reaction data class as a storage container for
 reading in the reaction mechanism being used and provides constructors
 for performing calculations with the data.  

 omega() returns the time rate of change of species c_num 
 
***************************************************************************
***************************************************************************/
// Reaction system flags
enum Reactions {NO_REACTIONS,
		CH4_1STEP,        // CH4 1 step HARDCODED mechanisms
		CH4_2STEP,        // CH4 2 step HARDCODED mechanisms
		C3H8_1STEP,       // C3H8 1 step HARDCODED mechanisms
		H2O2_1STEP,       // H2-O2 1 step HARDCODED mechanisms
		H2O2_2STEP,       // H2-O2 2 step HARDCODED mechanisms
		H2O2_8STEP,       // H2-O2 8 step HARDCODED mechanisms
		CH4_15STEP_ARM2,  // CH4 HARDCODED mechanisms based on GRI 2.11
		CH4_15STEP_ARM3,  // CH4 HARDCODED mechanisms based on GRI 3
		USER,             // User defined flag
		CANTERA};         // cantera flag

class Reaction_set{

private:  

  //Temp vectors used during omega & dSwdU
  double *kf; 
  double *kb; 
  double *M; 
  double *c; 
  double *c_denom; 
  double *r; 
  cplx   *cc;
  cplx   *rc;

protected:
public: 
  int reactset_flag;       //Reaction Set Flag
  React_data *reactions;   //each reaction rate data
  int num_reactions;       //number of reactions
  int num_species;         //number of total species
  int num_react_species;   //number of reacting species
  string *species;         //species used in reactions
  string Reaction_system;  //Reaction system name

  // cantera related objects
  string ct_mech_name;     //Reaction mechanism file path
  string ct_mech_file;     //Reaction mechanism file path
#ifdef _CANTERA_VERSION
  IdealGasMix* ct_gas;     //the Cantera IdealGasMix object
#endif

  // indexes of specific species
  static int iCO, iCO2, iH2O, iO2;

  Reaction_set(){ 
    reactset_flag=0; num_reactions=0; num_species=0; 
    num_react_species=0; reactions = NULL; species = NULL; 
    kf=NULL; kb=NULL; M=NULL; c=NULL; c_denom=NULL; r=NULL;
    cc=NULL; rc=NULL;
#ifdef _CANTERA_VERSION
    ct_gas=NULL;
#endif
  }
                      
  /******** Constructors *******************/
  //for hardcoded reactions
  void set_reactions(string &);
  //for user defined
  void set_species(string *, int);
  void set_reactions(int &,string*,double*,double*,double*);

  // cantera member functions
  void ct_load_mechanism(string &, string &);
  void ct_parse_mass_string( const string&, double* );
  void ct_parse_schmidt_string( const string&, double*);

  //setup storage after num_reactions & num_species set.
  void set_storage(void){
    kf = new double[num_reactions];        
    kb = new double[num_reactions];
    M  = new double[num_react_species];
    c  = new double[num_react_species];
    c_denom = new double[num_react_species]; 
    r  = new double[num_react_species]; 
    cc = new cplx[num_react_species];
    rc = new cplx[num_react_species]; 
  }

  //Operator Overloading 
  Reaction_set& operator =(const Reaction_set &W);

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const Reaction_set &W);
  //friend istream& operator >> (istream &in_file,  Reaction_set &W);

  // time rate change of the species concentration 
  template<class SOLN_pSTATE, class SOLN_cSTATE>
  void omega( SOLN_cSTATE &U, const SOLN_pSTATE &W, const int Flow_Type ) const;

  //Jacobian ( flag true for cfl calc)
  template<class SOLN_pSTATE, class SOLN_cSTATE>
  void dSwdU(DenseMatrix &dSwdU,const SOLN_pSTATE &W, const bool &CFL_flag, 
	     const int Flow_Type, const int Solver_type) const;

  //Jacobian computed using complex step method
  template<class SOLN_pSTATE>
  void Complex_Step_dSwdU(DenseMatrix &dSwdU, const SOLN_pSTATE &W) const;

  // jacobian 
  template<class SOLN_pSTATE, class SOLN_cSTATE>
  void Finite_Difference_dSwdU(DenseMatrix &dSwdU, const SOLN_pSTATE &W,
			       const int Flow_Type) const;

  void Deallocate();

  //destructor
  ~Reaction_set(){Deallocate();};
 
};


/**************** Destructor *******************************************/
inline void Reaction_set::Deallocate(){
  //deallocate memory
  if(reactions != NULL){  delete[] reactions; reactions = NULL;  }
  if(species != NULL){    delete[] species;   species = NULL;  }
  if(kf != NULL){  delete[] kf; kf = NULL;}
  if(kb != NULL){  delete[] kb; kb = NULL;}
  if(M != NULL){  delete[] M; M = NULL;}
  if(c != NULL){  delete[] c; c = NULL;}
  if(c_denom != NULL){  delete[] c_denom; c_denom = NULL;}
  if(r != NULL){  delete[] r; r = NULL;}
  if(cc!= NULL){ delete[] cc; cc = NULL;}
  if(rc!= NULL){ delete[] rc; rc = NULL;}
#ifdef _CANTERA_VERSION
  if(ct_gas != NULL){ delete ct_gas; ct_gas = NULL; }
#endif // _CANTERA_VERSION
}


/***************** Assignment ****************************************/
inline Reaction_set& Reaction_set::operator =(const Reaction_set &W){
  //self assignment protection
  if( this != &W){   

    // for all hard-coded and user cases
    if (W.reactset_flag != CANTERA) {
      string temp = W.Reaction_system;
      //copy assignment
      set_reactions(temp);
    }
    
    // for cantera case
    else {
      string mechanism_file = W.ct_mech_file;
      string mechanism_name = W.ct_mech_name;
      ct_load_mechanism(mechanism_file, mechanism_name);
    }

  } /* endif */
  return (*this);
}



/**************** I/O Operators ***************************************/
inline ostream &operator << (ostream &out_file, const Reaction_set &W) {
  out_file.setf(ios::scientific);
  out_file <<"\n "<<W.Reaction_system<<" "<<W.num_reactions
	   <<"  "<<W.num_species<<" "<<W.num_react_species<<endl;  
  //species
  for(int i=0; i<W.num_species; i++){
    out_file <<" "<<W.species[i];
  }
  //each reaction systems data
  if(W.reactset_flag != NO_REACTIONS){ 
    for(int i=0; i<W.num_reactions; i++){
      out_file <<"\n "<<W.reactions[i];
    }
  }
  out_file.unsetf(ios::scientific);
  return (out_file);
}

// istream &operator >> (istream &in_file, Reaction_set &W) {
//   in_file.setf(ios::skipws);
//   in_file >> W.react >> W.A >> W.b >> W.E;
//   in_file.unsetf(ios::skipws);
//   return (in_file);
// }


/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using the general law of mass action.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations
              of Chem2D_pState. (c_i*rho/M_i)   mol/m^3

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
inline void Reaction_set::omega(SOLN_cSTATE &U, const SOLN_pSTATE &W,  const int Flow_Type ) const{
 
  double Temp = W.T();  //K
  double Press= W.p;    // [Pa]
  double PressDyne = W.p*TEN; // N/m^2 -> dyne/cm^2
  double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3 
  double a,b, ans(ZERO);

  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND;  //kg/mol -> g/mol
    c[i] = W.spec[i].c;                        //unitless
  }

  switch(reactset_flag){
  case NO_REACTIONS:
    cerr<<"\n You shouldn't get here, NO_REACTIONS in Reaction_set::omeag(..)";
    exit(1);
    break;
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
  case CH4_1STEP: 
    //laminar case     
    a=1.0;//   a = 0.2; 
    b=1.0;//     b = 1.3;
    kf[0] = reactions[0].kf(Temp)*pow((W.SpecCon(0))/MILLION,a)*pow((W.SpecCon(1)/MILLION),b);        
    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //CH4
	ans = - kf[0];
	break;
      case 1 : //O2
	ans = - TWO*kf[0];
	break;
      case 2 : //CO2
	ans = kf[0];
	break;
      case 3 : //H2O
	ans = TWO*kf[0];
	break;
      };
      //ans in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3      
      U.rhospec[index].c = M[index]*ans*THOUSAND;  
    }
    break;
    
    //TWO STEP CH4
  case CH4_2STEP:
    kf[0] = reactions[0].kf(Temp)*pow(W.SpecCon(0)/MILLION,0.2)*pow(W.SpecCon(1)/MILLION,1.3);
    kb[0] = ZERO;
    kf[1] = reactions[1].kf(Temp)*(W.SpecCon(4)/MILLION)*pow(W.SpecCon(3)/MILLION,0.5)*pow(W.SpecCon(1)/MILLION,0.25); 
    kb[1] = reactions[1].kb(Temp)*(W.SpecCon(2)/MILLION);
   
    //cout<<"\n Sw "<<Temp<<" "<<reactions[0].kf(Temp)<<" "<<reactions[1].kf(Temp)<<" "<<reactions[1].kb(Temp);

    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //CH4
	ans = - kf[0];
	break;
      case 1:  //O2
	ans = -1.5*kf[0] - HALF*kf[1] + HALF*kb[1];
	break;
      case 2:  //CO2
	ans = kf[1]  - kb[1];
	break;
      case 3:  //H2O
	ans = TWO*kf[0];
	break;
      case 4:  //CO
	ans = kf[0] - kf[1] + kb[1];
 	break;
      };     
      //ans in kg/m^3*s     g/mol *(mol/cm^3*s)*1e3       
      U.rhospec[index].c = M[index]*ans*THOUSAND;
    }
    break;
  
  // ONE STEP C3H8 mechanism
  case C3H8_1STEP:
    cerr<<"\nC3H8_1STEP not implemented.";
    exit(1);
    break;


    //ONE STEP H2&O2
  case H2O2_1STEP:
   
      kf[0] = reactions[0].kf(Temp);
      
      for(int index =0; index<num_react_species; index++){
	switch(index) {
	case 0 : //H2
	  ans = - TWO*kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	case 1 : //O2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	case 2 : //H2O
	  ans = TWO*kf[0]*W.SpecCon(0)*W.SpecCon(0)*W.SpecCon(1);
	  break;
	};
	//ans in kg/m^3*s   
	U.rhospec[index].c = M[index]*ans*THOUSAND;
      }
      break;
  
      //TWO STEP H2&O2
  case H2O2_2STEP:
   
      kf[0] = reactions[0].kf(Temp,W.spec[0].c,W.spec[1].c,W.spec[4].c);
      kb[0] = kf[0]/(reactions[0].keq(W,Temp));
      kf[1] = reactions[1].kf(Temp,W.spec[0].c,W.spec[1].c,W.spec[4].c);
      kb[1] = kf[1]/(reactions[1].keq(W,Temp));
      //cout<<"\n kf1 "<< kf[0]<<" kb1 "<< kb[0]<<" kf2 "<<kf[1]<<" kb2 "<< kb[1];
      
      for(int index =0; index<num_react_species; index++){
	switch(index) {
	case 0 : //H2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 + kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12
	    - kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 + kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12;
	  break;
	case 1 : //O2
	  ans = - kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 + kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12;
	  break;
	case 2 : //OH
	  ans = TWO*(kf[0]*W.SpecCon(0)*W.SpecCon(1)*1e-12 - kb[0]*W.SpecCon(2)*W.SpecCon(2)*1e-12)
	    + TWO*( -kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 + kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12);
	  break;
	case 3 : //H2O
	  ans = TWO*(kf[1]*W.SpecCon(2)*W.SpecCon(2)*W.SpecCon(0)*1e-18 - kb[1]*W.SpecCon(3)*W.SpecCon(3)*1e-12);      
	  break;
	};
	//ans in kg/m^3*s   
	U.rhospec[index].c = M[index]*ans*THOUSAND;
      }
      break;

    // 8 STEP H2 & O2
  case H2O2_8STEP: 
    
    for(int spec =0; spec<num_reactions; spec++){
      kf[spec] = reactions[spec].kf(Temp);
      kb[spec] = kf[spec]/reactions[spec].keq(W,Temp);
    }
 
    for(int index =0; index<num_react_species; index++){
      switch(index) {
      case 0 : //O
	ans = (kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]-kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]-kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]+kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]+kf[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])-kb[3]*rho*rho*c[4]/M[4]*c[0]/M[0]-kf[6]
			*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]+kb[6]*rho*rho*c[5]/M[5]*c[4]/M[4]
			-2.0*kf[7]*rho*rho*rho*c[0]*c[0]/(M[0]*M[0])*c[4]/M[4]+2.0*kb[7]*rho*rho*c[1]/M
			[1]*c[4]/M[4]);
	
	break;
      case 1 : //O2
	ans = (-kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]+kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[7]*rho*rho*rho*c[0]*c[0]/(M[0]*M[0])*c[4]/M[4]-kb[7]*rho*rho*c
			[1]/M[1]*c[4]/M[4]);
  
	break;
      case 2 : //H
	ans = (-kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]+kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]-kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]+kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]-kb[2]*rho*rho*c[2]/M[2]*c[4]/M[4]-2.0*kf
			[4]*rho*rho*rho*c[2]*c[2]/(M[2]*M[2])*c[4]/M[4]+2.0*kb[4]*rho*rho*c[3]/M[3]*c
			[4]/M[4]-kf[5]*rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]+kb[5]*rho*rho*c[4]*c
			 [4]/(M[4]*M[4])-kf[6]*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]+kb[6]*rho*rho*c
			[5]/M[5]*c[4]/M[4]);	
	break;
      case 3 : //H2
	ans = (-kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]+kb[1]*rho*rho*c[5]/M
			 [5]*c[2]/M[2]-kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]+kb[2]*rho*rho*c[2]/M[2]*c[4]/M
			 [4]+kf[4]*rho*rho*rho*c[2]*c[2]/(M[2]*M[2])*c[4]/M[4]-kb[4]*rho*rho*c[3]/M[3]*c
			 [4]/M[4]);
	break;
	
      case 4 : //H2O
	ans = (kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]-kb[2]*rho*rho*c[2]/M
			[2]*c[4]/M[4]+kf[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])-kb[3]*rho*rho*c[4]/M[4]*c[0]/
			M[0]+kf[5]*rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]-kb[5]*rho*rho*c[4]*c[4]/(M[4]*M[4]));
	break;

      case 5 : // OH
	ans = (kf[0]*rho*rho*c[2]/M[2]*c[1]/M[1]-kb[0]*rho*rho*c[5]/M
			[5]*c[0]/M[0]+kf[1]*rho*rho*c[0]/M[0]*c[3]/M[3]-kb[1]*rho*rho*c[5]/M[5]*c[2]/M
			[2]-kf[2]*rho*rho*c[5]/M[5]*c[3]/M[3]+kb[2]*rho*rho*c[2]/M[2]*c[4]/M[4]-2.0*kf
			[3]*rho*rho*c[5]*c[5]/(M[5]*M[5])+2.0*kb[3]*rho*rho*c[4]/M[4]*c[0]/M[0]-kf[5]*
			rho*rho*rho*c[2]/M[2]*c[5]/M[5]*c[4]/M[4]+kb[5]*rho*rho*c[4]*c[4]/(M[4]*M[4])+
			kf[6]*rho*rho*rho*c[0]/M[0]*c[2]/M[2]*c[4]/M[4]-kb[6]*rho*rho*c[5]/M[5]*c[4]/M
			[4]);
	break;
      };
      U.rhospec[index].c = M[index]*ans*THOUSAND;
   }
   break;


    // 15 step CH4 based on GRI 2.11
  case CH4_15STEP_ARM2:
    // r(Wdot) - mol/(cm^3*s)

    // compute reaction rates calling the subroutine CKWYP15STEP211
    ckwyp15step211_(PressDyne, Temp, c, r);
    for (int i=0; i<num_react_species; ++i) {
      // in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3             
      U.rhospec[i].c = M[i]*r[i]*THOUSAND;
    }
    break;


    // 15 step CH4 based on GRI 3
  case CH4_15STEP_ARM3:
    // r(Wdot) - mol/(cm^3*s)

    // compute reaction rates calling the subroutine CKWYP15STEP3
    ckwyp15step30_(PressDyne, Temp, c, r);
    for (int i=0; i<num_react_species; ++i) {
      // in kg/m^3*s   g/mol *(mol/cm^3*s)*1e3             
      U.rhospec[i].c = M[i]*r[i]*THOUSAND;
    }
    break;


  //---------------------------------//
  //------------ CANTERA ------------//
  //---------------------------------//
#ifdef _CANTERA_VERSION

  case CANTERA:

    // set the gas state
    ct_gas->setState_TPY(Temp, Press, c);

    // get the species net production rates
    ct_gas->getNetProductionRates(r);

    // set the Chem2D_cState production rates
    for (int i=0; i<num_react_species; i++){
      U.rhospec[i].c = r[i]*M[i];
    }

    break;

#endif // _CANTERA_VERSION

    //---------------------------------//
    //----- User Specified ------------//
    //---------------------------------//
  case USER:
    cerr<<"\nUser specified not currently available in shareware version";
    exit(1);
    break;
  default:
    break;
  };


} //end omega()

/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables.  Ie it returns conserved data.
  
   dSwdU:  Matrix of source terms 
   W:      Primitive State data.

  These are added to the RHS so += is used.

  NOTE: Be carful of mass fractions in the denominator as 
        if they go to ZERO, it will cause floating point exceptions
       (ie. Divison by zero);  
************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
inline void Reaction_set::dSwdU(DenseMatrix &dSwdU, const SOLN_pSTATE &W, 
				const bool &CFL_flag, const int Flow_Type,
				const int solver_type) const{

  /***************** Local Variables *****************/
  double Temp = W.T();
  double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3


//   double Con0, Con1, Rt;
//   double dkf0_dp, dkf0_drho;
//   double dkf0_dc1, dkf0_dc2, dkf0_dc3, dkf0_dc4; 
//   //BEING LAZY INSTEAD OF ANALYTICALY DETERMINING dSdU
//   DenseMatrix dSwdW(W.NUM_VAR_CHEM2D-1,W.NUM_VAR_CHEM2D-1,ZERO);  //SHOULD BE MOVED TO TEMP AS WELL!!!
//   DenseMatrix dWdQ(W.NUM_VAR_CHEM2D-1,W.NUM_VAR_CHEM2D-1,ZERO);
  
  double VALUE = TOLER; //sqrt(TOLER);
  //////////////////////////////////////////////////////////////////////////////
  // THIS HACK DOESN'T REALLY WORK, ESPECIALLY FOR 2STEP EQUATIONS !!!!!!!!!!!!!
  //////////////////////////////////////////////////////////////////////////////
  if(solver_type == IMPLICIT){ VALUE=TOLER*TOLER; }

  // intitialize
  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND; //kg/mol -> g/mol 
    c[i] = W.spec[i].c;
    
    //For handling ~= ZERO mass fractions that appear in the denominator of dSwdU
    //by setting a lower tolerance allowed, and if below that set to tolerance 
    if( c[i] < VALUE){
      c_denom[i] = VALUE;
    } else {
      c_denom[i] = c[i];
    }
  }

  int NUM_VAR = W.NumVarSansSpecies();

  /*******************************************
   *  Reaction Mechanism Jacobians           *
   *                                         * 
   *******************************************/
  switch(reactset_flag){
  case NO_REACTIONS:
    //this case shouldn't be called
    cerr<<" Ummm, this dSwdU shouldn't be called for with NO_REACTIONS, check quad_singleblock.";
    break;
    
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
    // This is far from an elegant solution, but its easily created using maple.
    // It can probably be simplified for faster computation.

  case CH4_1STEP:
    kf[0] = reactions[0].kf(Temp); 

    // One forward reaction ... so only depending on the fuel and oxidize ..  
    // for coef a=0.2, b=1.3
 //    dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1);
//     dSwdU(NUM_VAR+1,NUM_VAR) += -0.4*M[1]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
//     dSwdU(NUM_VAR+2,NUM_VAR) += 0.2*M[2]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
//     dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];    
//     dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
//     dSwdU(NUM_VAR+1,NUM_VAR+1) += -0.26E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3);
//     dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.13E1*M[2]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
//     dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];  
 
    //for coef a=b=1.0
    dSwdU(NUM_VAR,NUM_VAR) += -1.0*kf[0]*rho*c[1]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+1) += -1.0*kf[0]*rho*c[0]/M[1];
    dSwdU(NUM_VAR+1,NUM_VAR) += -2.0*kf[0]*rho*c[1]/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -2.0*kf[0]*rho*c[0]/M[0];    
    dSwdU(NUM_VAR+2,NUM_VAR) +=  M[2]*kf[0]*rho*c[1]/M[1]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += M[2]*kf[0]*rho*c[0]/M[0]/M[1];
    dSwdU(NUM_VAR+3,NUM_VAR) += 2.0*M[3]*kf[0]*rho*c[1]/M[1]/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+1) += 2.0*M[3]*kf[0]*rho*c[0]/M[0]/M[1];


//     /********** dSwdW ***********************/
//     Con0 = rho*c[0]/M[0];
//     Con1 = rho*c[1]/M[0];
    
//     Rt = W.Rtot(); 
//     dkf0_dp = reactions[0].dkf_dT(Temp)/(rho*Rt); ///THOUSAND;
//     dkf0_drho = reactions[0].dkf_dT(Temp)*(-Temp/rho);
    
//     //take into account kf(rho,p,and ci)
//     dkf0_dc1 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[0].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc2 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[1].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc3 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[2].Rs()-W.specdata[4].Rs()))*Con0*Con1;
//     dkf0_dc4 =  reactions[0].dkf_dT(Temp)*(-Temp/Rt*(W.specdata[3].Rs()-W.specdata[4].Rs()))*Con0*Con1;

//     dSwdW(NUM_VAR,NUM_VAR)  += -M[0]*dkf0_dc1;
//     dSwdW(NUM_VAR+1,NUM_VAR) += -TWO*M[1]*dkf0_dc1;
//     dSwdW(NUM_VAR+2,NUM_VAR) += M[2]*dkf0_dc1;
//     dSwdW(NUM_VAR+3,NUM_VAR) += TWO*M[3]*dkf0_dc1;

//     dSwdW(NUM_VAR,NUM_VAR+1)  += -M[0]*dkf0_dc2;
//     dSwdW(NUM_VAR+1,NUM_VAR+1) += -TWO*M[1]*dkf0_dc2;
//     dSwdW(NUM_VAR+2,NUM_VAR+1) += M[2]*dkf0_dc2;
//     dSwdW(NUM_VAR+3,NUM_VAR+1) += TWO*M[3]*dkf0_dc2;
    
//     dSwdW(NUM_VAR,NUM_VAR+2)  += -M[0]*dkf0_dc3;
//     dSwdW(NUM_VAR+1,NUM_VAR+2) += -TWO*M[1]*dkf0_dc3;
//     dSwdW(NUM_VAR+2,NUM_VAR+2) += M[2]*dkf0_dc3;
//     dSwdW(NUM_VAR+3,NUM_VAR+2) += TWO*M[3]*dkf0_dc3;

//     dSwdW(NUM_VAR,NUM_VAR+3)  += -M[0]*dkf0_dc4;
//     dSwdW(NUM_VAR+1,NUM_VAR+3) += -TWO*M[1]*dkf0_dc4;
//     dSwdW(NUM_VAR+2,NUM_VAR+3) += M[2]*dkf0_dc4;
//     dSwdW(NUM_VAR+3,NUM_VAR+3) += TWO*M[3]*dkf0_dc4;

//     dSwdW(NUM_VAR,0) += -M[0]*dkf0_drho*Con0*Con1 - kf[0]*Con1*c[0] - M[0]*kf[0]*Con0*c[1]/M[1];    
//     dSwdW(NUM_VAR+1,0) += -TWO*M[1]*dkf0_drho*Con0*Con1 - TWO*kf[0]*M[1]*Con1*c[0]/M[0] - TWO*kf[0]*Con0*c[1];    
//     dSwdW(NUM_VAR+2,0) += M[2]*dkf0_drho*Con0*Con1 + kf[0]*M[2]*Con1*c[0]/M[0] + M[2]*kf[0]*Con0*c[1]/M[1];    
//     dSwdW(NUM_VAR+3,0) += TWO*M[3]*dkf0_drho*Con0*Con1 + TWO*kf[0]*M[3]*Con1*c[0]/M[0] + TWO*M[3]*kf[0]*Con0*c[1]/M[1];   

//     dSwdW(NUM_VAR,3) += -M[0]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+1,3) += -TWO*M[1]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+2,3) += M[2]*dkf0_dp*Con0*Con1;
//     dSwdW(NUM_VAR+3,3) += TWO*M[3]*dkf0_dp*Con0*Con1;
   
//     dSwdW(NUM_VAR,NUM_VAR)   += -kf[0]*Con1*rho*THOUSAND;
//     dSwdW(NUM_VAR+1,NUM_VAR) += -TWO*M[1]*kf[0]*Con1*rho/M[0]*THOUSAND;
//     dSwdW(NUM_VAR+2,NUM_VAR) += M[2]*kf[0]*Con1*rho/M[0]*THOUSAND;
//     dSwdW(NUM_VAR+3,NUM_VAR) += TWO*M[3]*kf[0]*Con1*rho/M[0]*THOUSAND;

//     dSwdW(NUM_VAR,NUM_VAR+1)   += -M[0]*kf[0]*Con0*rho/M[1]*THOUSAND;
//     dSwdW(NUM_VAR+1,NUM_VAR+1) += -TWO*kf[0]*Con0*rho*THOUSAND;
//     dSwdW(NUM_VAR+2,NUM_VAR+1) += M[2]*kf[0]*Con0*rho/M[1]*THOUSAND;
//     dSwdW(NUM_VAR+3,NUM_VAR+1) += TWO*M[3]*kf[0]*Con0*rho/M[1]*THOUSAND;

//     // LAZY dSwdU = dSwdW * dWdU
//     W.dWdU(dWdQ,Flow_Type); 
//     dSwdU += dSwdW*dWdQ;
//     /***************************************/

    //this is a work around for the delta t calculation using an unesseccarily small value
    //when c[0] -> ZERO
    if(c_denom[0] != c[0] && CFL_flag){ dSwdU(NUM_VAR,NUM_VAR)=ZERO; }
    
    break;


    //TWO STEP CH4   
  case CH4_2STEP:  
    /******************** ORIGINAL ****************************************/
    //still some issues with units ??? ie.  dSwdU(6,6) ??
    //which is also on the diagonal so messes with CFL???
    kf[0] = reactions[0].kf(Temp);      
    kb[0] = ZERO;
    kf[1] = reactions[1].kf(Temp);  
    kb[1] = reactions[1].kb(Temp);

    //cout<<"\n dSwdU "<<Temp<<" "<<reactions[0].kf(Temp)<<" "<<reactions[1].kf(Temp)<<" "<<reactions[1].kb(Temp);

    dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3);
    dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];      
    
    dSwdU(NUM_VAR+1,NUM_VAR) += -0.3*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -rho*(0.195E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
			*M[4]+0.125*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/M[1]/M[4]/pow(rho*c_denom[1]/M[1],0.75);
    dSwdU(NUM_VAR+1,NUM_VAR+2) += 0.5*M[1]*kb[1]/M[2];
    dSwdU(NUM_VAR+1,NUM_VAR+3) += -0.25*M[1]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+1,NUM_VAR+4) += -0.5*M[1]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
   
    dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.25*M[2]*kf[1]*rho*c[4]/M[4]*sqrt(rho*c[3]/M[3])/pow(rho*c_denom[1]/M[1],0.75)/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -kb[1];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += 0.5*M[2]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+2,NUM_VAR+4) += M[2]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
    
    dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];    
   
    dSwdU(NUM_VAR+4,NUM_VAR) += 0.2*M[4]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
    dSwdU(NUM_VAR+4,NUM_VAR+1) += rho*(0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
		       *M[4]-0.25*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/(M[1]*M[1])/pow(rho*c_denom[1]/M[1],0.75);
    dSwdU(NUM_VAR+4,NUM_VAR+2) += M[4]*kb[1]/M[2];
    dSwdU(NUM_VAR+4,NUM_VAR+3) += -0.5*kf[1]*rho*c[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
    dSwdU(NUM_VAR+4,NUM_VAR+4) += -1.0*kf[1]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);

//     if(!CFL_flag) cout<<"\n ORIGINAL\n"<<dSwdU;
//     dSwdU.zero();

//     /******************** MODIFIED ****************************************/   
//     //with kf(T), kb(T)

//     kf[0] = reactions[0].kf(Temp);      //UNITS ISSUES due to rho and M and kf's
//     kb[0] = ZERO;
//     kf[1] = reactions[1].kf(Temp);  
//     kb[1] = reactions[1].kb(Temp);
    
//     double dkf_dT[2], dkb_dT[2];
//     dkf_dT[0] = reactions[0].dkf_dT(Temp);
//     dkb_dT[0] = ZERO;
//     dkf_dT[1] = reactions[1].dkf_dT(Temp);
//     dkb_dT[1] = reactions[1].dkb_dT(Temp);

//     double T = Temp;
//     double Cp = W.Cp();      
//     double u = W.v.x*HUNDRED; 
//     double v = W.v.y*HUNDRED; 
//     //double p = W.p;      
//     double Rtot = W.Rtot();   
//     double htot = W.h();      
//     double h[6],Rs[6];

//     for(int i=0; i<num_species; i++){    
//       h[i] = (W.specdata[i].Enthalpy(Temp)+W.specdata[i].Heatofform());  
//       Rs[i] = (W.specdata[i].Rs()); 
//     }    

//       if(!CFL_flag) cout<<"\n NEW\n"<<dSwdU;
   
    //this is a work around for the delta t calculation using an unesseccarily small value
    //when CH4 & O2 -> ZERO
    //if( (c[0] < TOLER && CFL_flag) || (c[0] < TOLER && CFL_flag)  ){  
    
    if( c_denom[0] != c[0] && CFL_flag ){
      for(int i=0; i<num_react_species; i++){
	dSwdU(NUM_VAR+i,NUM_VAR+i)=ZERO; 
      }
    }
    break;
  
    //take the minum value between laminar and turbulent concentrations
    //   if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
    // 	Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){
    //     }


  // ONE STEP C3H8 mechanism
  case C3H8_1STEP:
    cerr<<"\nC3H8_1STEP not implemented.";
    exit(1);
    break;


  case H2O2_1STEP:  
    kf[0] = reactions[0].kf(Temp);
    rho = rho;
    M[0] = M[0];
    M[1] = M[1];

    dSwdU(NUM_VAR,NUM_VAR) += -4.0*rho*rho/M[0]*kf[0]*c[0]*c[1]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+1) = -2.0*rho*rho/M[0]*kf[0]*c[0]*c[0]/M[1];
    
    dSwdU(NUM_VAR+1,NUM_VAR) = -2.0*rho*rho*kf[0]*c[0]/(M[0]*M[0])*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+1) = -rho*rho*kf[0]*c[0]*c[0]/(M[0]*M[0]);
    
    dSwdU(NUM_VAR+2,NUM_VAR) = 4.0*rho*rho*M[2]*kf[0]*c[0]/(M[0]*M[0])*c[1]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+1)= 2.0*rho*rho*M[2]*kf[0]*c[0]*c[0]/(M[0]*M[0])/M[1];
  
    break;
   //TWO STEP H2 and O2
  case H2O2_2STEP:
    kf[0] = reactions[0].kf(Temp,c[0],c[1],W.spec[4].c);
    kb[0] = kf[0]/reactions[0].keq(W,Temp);
    kf[1] = reactions[1].kf(Temp,c[0],c[1],W.spec[4].c);
    kb[1] = kf[1]/reactions[1].keq(W,Temp);  
  
    dSwdU(NUM_VAR,NUM_VAR) += -rho*(kf[0]*c[1]*M[2]*M[2]+kf[1]*rho*c[2]*c[2]*M[1])/M[1]/(M[2]*M[2]);
    dSwdU(NUM_VAR,NUM_VAR+1) += -kf[0]*rho*c[0]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+2) += -2.0*rho*c[2]*(-kb[0]*M[0]+kf[1]*rho*c[0])/(M[2]*M[2]);
    dSwdU(NUM_VAR,NUM_VAR+3) += 2.0*M[0]*kb[1]*rho*c[3]/(M[3]*M[3]);    
  
    dSwdU(NUM_VAR+1,NUM_VAR) += -kf[0]*rho/M[0]*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -kf[0]*rho*c[0]/M[0];
    dSwdU(NUM_VAR+1,NUM_VAR+2) += 2.0*M[1]*kb[0]*rho*c[2]/(M[2]*M[2]);    
   
    dSwdU(NUM_VAR+2,NUM_VAR) += rho*(0.2E1*kf[0]*c[1]*M[2]*M[2]-0.2E1*kf[1]*rho*c[2]*c[2]*M[1])/M[2]/M[0]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += 2.0*M[2]*kf[0]*rho*c[0]/M[0]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -rho*c[2]*(0.4E1*kb[0]*M[0]+0.4E1*kf[1]*rho*c[0])/M[2]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += 4.0*M[2]*kb[1]*rho*c[3]/(M[3]*M[3]);    
   
    dSwdU(NUM_VAR+3,NUM_VAR) += 2.0*M[3]*kf[1]*rho*rho*c[2]*c[2]/(M[2]*M[2])/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+2) += 4.0*M[3]*kf[1]*rho*rho*c[2]/(M[2]*M[2])*c[0]/M[0];
    dSwdU(NUM_VAR+3,NUM_VAR+3) += -4.0/M[3]*kb[1]*rho*c[3];  

    break;
    // 8 STEP H2 & O2
  case H2O2_8STEP:  

    for(int spec =0; spec<num_reactions; spec++){
      kf[spec] = reactions[spec].kf(Temp);
      kb[spec] = kf[spec]/reactions[spec].keq(W,Temp); 
    }

    dSwdU(NUM_VAR,NUM_VAR) += -1/M[0]*rho*(kb[0]*c[5]*M[0]*M[3]*M[4]*M[2]+kf[1]*c[3]*M[5]
					   *M[0]*M[4]*M[2]+kb[3]*c[4]*M[5]*M[0]*M[3]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[0]*M
			       [3]+4.0*kf[7]*rho*c[0]*c[4]*M[5]*M[3]*M[2])/M[5]/M[3]/M[4]/M[2];
    dSwdU(NUM_VAR,NUM_VAR+1) += M[0]*rho*(kf[0]*c[2]*M[4]+2.0*kb[7]*c[4]*M[2])/M[2]/M[1]/M[4];
    dSwdU(NUM_VAR,NUM_VAR+2) += -rho*(-kf[0]*c[1]*M[5]*M[0]*M[4]-kb[1]*c[5]*M[1]*M[0]*M[4]+
			kf[6]*rho*c[0]*c[4]*M[1]*M[5])/M[2]/M[1]/M[5]/M[4];
    dSwdU(NUM_VAR,NUM_VAR+3) += -kf[1]*rho*c[0]/M[3];
    dSwdU(NUM_VAR,NUM_VAR+4) += -1/M[0]*rho*(kb[3]*c[0]*M[0]*M[2]*M[5]*M[1]+kf[6]*rho*c[0]*
			       c[2]*M[0]*M[5]*M[1]-kb[6]*c[5]*M[0]*M[0]*M[2]*M[1]+2.0*kf[7]*rho*c[0]*c[0]*M[2]
			       *M[5]*M[1]-2.0*kb[7]*c[1]*M[0]*M[0]*M[2]*M[5])/M[4]/M[2]/M[5]/M[1];
    dSwdU(NUM_VAR,NUM_VAR+5) += rho*(-kb[0]*c[0]*M[5]*M[2]*M[4]+kb[1]*c[2]*M[5]*M[0]*M[4]+
		       2.0*kf[3]*c[5]*M[0]*M[2]*M[4]+kb[6]*c[4]*M[5]*M[0]*M[2])/(M[5]*M[5])/M[2]/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR) += M[1]*rho*(kb[0]*c[5]*M[0]*M[4]+2.0*kf[7]*rho*c[0]*c[4]*M[5])/M[5]/(M[0]*M[0])/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+1) += -rho*(kf[0]*c[2]*M[4]+kb[7]*c[4]*M[2])/M[2]/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+2) += -kf[0]*rho/M[2]*c[1];
    dSwdU(NUM_VAR+1,NUM_VAR+4) += rho*(kf[7]*rho*c[0]*c[0]*M[1]-kb[7]*c[1]*M[0]*M[0])/(M[0]*M[0])/M[4];
    dSwdU(NUM_VAR+1,NUM_VAR+5) += M[1]*kb[0]*rho/M[5]*c[0]/M[0];
   
    dSwdU(NUM_VAR+2,NUM_VAR) += -rho*(-kb[0]*c[5]*M[3]*M[4]*M[2]-kf[1]*c[3]*M[5]*M[4]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[3])/M[0]/M[5]/M[3]/M[4];
    dSwdU(NUM_VAR+2,NUM_VAR+1) += -kf[0]*rho*c[2]/M[1];
    dSwdU(NUM_VAR+2,NUM_VAR+2) += -1/M[2]*rho*(kf[0]*c[1]*M[2]*M[5]*M[4]*M[0]+kb[1]*c[5]*M[2]
			       *M[1]*M[4]*M[0]+kb[2]*c[4]*M[2]*M[1]*M[5]*M[0]+4.0*kf[4]*rho*c[2]*c[4]*M[1]*M
			       [5]*M[0]+kf[5]*rho*c[5]*c[4]*M[2]*M[1]*M[0]+kf[6]*rho*c[0]*c[4]*M[2]*M[1]*M[5])/M[1]/M[5]/M[4]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+3) += M[2]*rho*(kf[1]*c[0]*M[5]*M[4]+kf[2]*c[5]*M[0]*M[4]+2.0*kb
			    [4]*c[4]*M[0]*M[5])/M[0]/M[3]/M[5]/M[4];
    dSwdU(NUM_VAR+2,NUM_VAR+4) += -1/M[2]*rho*(kb[2]*c[2]*M[2]*M[4]*M[3]*M[5]*M[0]+2.0*kf[4]*
			       rho*c[2]*c[2]*M[4]*M[3]*M[5]*M[0]-2.0*kb[4]*c[3]*M[2]*M[2]*M[4]*M[5]*M[0]+kf[5]
			       *rho*c[2]*c[5]*M[2]*M[4]*M[3]*M[0]-2.0*kb[5]*c[4]*M[2]*M[2]*M[3]*M[5]*M[0]+kf
			       [6]*rho*c[0]*c[2]*M[2]*M[4]*M[3]*M[5]-kb[6]*c[5]*M[2]*M[2]*M[4]*M[3]*M[0])/(M[4]*M[4])/M[3]/M[5]/M[0];
    dSwdU(NUM_VAR+2,NUM_VAR+5) += -rho*(-kb[0]*c[0]*M[3]*M[4]*M[2]+kb[1]*c[2]*M[0]*M[3]*M[4]-
			kf[2]*c[3]*M[0]*M[4]*M[2]+kf[5]*rho*c[2]*c[4]*M[0]*M[3]-kb[6]*c[4]*M[0]*M[3]*M
			  [2])/M[0]/M[5]/M[3]/M[4]; 
    dSwdU(NUM_VAR+3,NUM_VAR) += -kf[1]*rho/M[0]*c[3];
    dSwdU(NUM_VAR+3,NUM_VAR+2) += M[3]*rho*(kb[1]*c[5]*M[2]*M[4]+kb[2]*c[4]*M[5]*M[2]+2.0*kf
			    [4]*rho*c[2]*c[4]*M[5])/M[5]/(M[2]*M[2])/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+3) += -rho*(kf[1]*c[0]*M[5]*M[4]+kf[2]*c[5]*M[0]*M[4]+kb[4]*c[4]*
			M[0]*M[5])/M[0]/M[5]/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+4) += rho*(kb[2]*c[2]*M[2]*M[3]+kf[4]*rho*c[2]*c[2]*M[3]-kb[4]*c
		       [3]*M[2]*M[2])/(M[2]*M[2])/M[4];
    dSwdU(NUM_VAR+3,NUM_VAR+5) += -rho*(-kb[1]*c[2]*M[3]+kf[2]*c[3]*M[2])/M[5]/M[2];
   
    dSwdU(NUM_VAR+4,NUM_VAR) += -kb[3]*rho*c[4]/M[0];
    dSwdU(NUM_VAR+4,NUM_VAR+2) += rho*c[4]*(-kb[2]*M[5]+kf[5]*rho*c[5])/M[2]/M[5];
    dSwdU(NUM_VAR+4,NUM_VAR+3) += M[4]*kf[2]*rho*c[5]/M[5]/M[3];
    dSwdU(NUM_VAR+4,NUM_VAR+4) += 1/M[4]*rho*(-kb[2]*c[2]*M[4]*M[0]*M[5]-kb[3]*c[0]*M[2]*M[4]
			      *M[5]+kf[5]*rho*c[2]*c[5]*M[4]*M[0]-2.0*kb[5]*c[4]*M[2]*M[0]*M[5])/M[2]/M[0]/M[5];
    dSwdU(NUM_VAR+4,NUM_VAR+5) += rho*(kf[2]*c[3]*M[5]*M[2]*M[4]+2.0*kf[3]*c[5]*M[3]*M[2]*M
		       [4]+kf[5]*rho*c[2]*c[4]*M[5]*M[3])/(M[5]*M[5])/M[3]/M[2];
    dSwdU(NUM_VAR+5,NUM_VAR) += rho*(-kb[0]*c[5]*M[3]*M[4]*M[2]+kf[1]*c[3]*M[5]*M[4]*M[2]+
		       2.0*kb[3]*c[4]*M[5]*M[3]*M[2]+kf[6]*rho*c[2]*c[4]*M[5]*M[3])/M[0]/M[3]/M[4]/M[2];
    dSwdU(NUM_VAR+5,NUM_VAR+1) += M[5]*kf[0]*rho*c[2]/M[2]/M[1];
    dSwdU(NUM_VAR+5,NUM_VAR+2) += -rho*(-kf[0]*c[1]*M[5]*M[0]*M[4]+kb[1]*c[5]*M[1]*M[0]*M[4]-
			kb[2]*c[4]*M[1]*M[5]*M[0]+kf[5]*rho*c[5]*c[4]*M[1]*M[0]-kf[6]*rho*c[0]*c[4]*M
			[1]*M[5])/M[2]/M[1]/M[0]/M[4];
    dSwdU(NUM_VAR+5,NUM_VAR+3) += -rho*(-kf[1]*c[0]*M[5]+kf[2]*c[5]*M[0])/M[0]/M[3];
    dSwdU(NUM_VAR+5,NUM_VAR+4) += -rho*(-kb[2]*c[2]*M[4]*M[0]*M[5]-2.0*kb[3]*c[0]*M[2]*M[4]*M
			[5]+kf[5]*rho*c[2]*c[5]*M[4]*M[0]-2.0*kb[5]*c[4]*M[2]*M[0]*M[5]-kf[6]*rho*c[0]*
			c[2]*M[4]*M[5]+kb[6]*c[5]*M[2]*M[4]*M[0])/M[2]/(M[4]*M[4])/M[0];
    dSwdU(NUM_VAR+5,NUM_VAR+5) += -1/M[5]*rho*(kb[0]*c[0]*M[5]*M[2]*M[3]*M[4]+kb[1]*c[2]*M[5]
			       *M[0]*M[3]*M[4]+kf[2]*c[3]*M[5]*M[0]*M[2]*M[4]+4.0*kf[3]*c[5]*M[0]*M[2]*M[3]*M
			       [4]+kf[5]*rho*c[2]*c[4]*M[5]*M[0]*M[3]+kb[6]*c[4]*M[5]*M[0]*M[2]*M[3])/M[0]/M[2]/M[3]/M[4];    
      
    break;


    // 15 step CH4 based on GRI 2.11
  case CH4_15STEP_ARM2:
    Complex_Step_dSwdU<SOLN_pSTATE>(dSwdU, W);
    
    //cout << "\n Complex Step dSwdU: " << dSwdU;
    //dSwdU.zero();
    //Finite_Difference_dSwdU(dSwdU, W, Flow_Type, Simple_Chemistry);
    //cout << "\n Finite Difference dSwdU: " << dSwdU;
    break;


    // 15 step CH4 based on GRI 3
  case CH4_15STEP_ARM3:
    Complex_Step_dSwdU<SOLN_pSTATE>(dSwdU, W);
    break;

  //---------------------------------//
  //------------ CANTERA ------------//
  //---------------------------------//
#ifdef _CANTERA_VERSION

  case CANTERA:
    Finite_Difference_dSwdU<SOLN_pSTATE,SOLN_cSTATE>(dSwdU, W, Flow_Type);
    break;

#endif // _CANTERA_VERSION

    //---------------------------------//
    //----- User Specified ------------//
    //---------------------------------//
  case USER:
    cerr<<"\nUser specified not set up yet";
    exit(1);
    break;
  default:
    //Do nothing (i.e. Jacobian = ZERO)
    break;
  };

    
} //end dSwdU



    /******************** ORIGINAL ****************************************/
  //TWO STEP CH4   
//  case CH4_2STEP:  
//     //still some issues with units ??? ie.  dSwdU(6,6) ??
//     //which is also on the diagonal so messes with CFL???
//     kf[0] = reactions[0].kf(Temp);      
//     kb[0] = ZERO;
//     kf[1] = reactions[1].kf(Temp);  
//     kb[1] = reactions[1].kb(Temp);

//     dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3);
//     dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];      
    
//     dSwdU(NUM_VAR+1,NUM_VAR) += -0.3*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]*pow(rho*c[1]/M[1],0.3)/M[0];
//     dSwdU(NUM_VAR+1,NUM_VAR+1) += -rho*(0.195E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
// 			*M[4]+0.125*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/M[1]/M[4]/pow(rho*c_denom[1]/M[1],0.75);
//     dSwdU(NUM_VAR+1,NUM_VAR+2) += 0.5*M[1]*kb[1]/M[2];
//     dSwdU(NUM_VAR+1,NUM_VAR+3) += -0.25*M[1]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
//     dSwdU(NUM_VAR+1,NUM_VAR+4) += -0.5*M[1]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
   
//     dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.25*M[2]*kf[1]*rho*c[4]/M[4]*sqrt(rho*c[3]/M[3])/pow(rho*c_denom[1]/M[1],0.75)/M[1];
//     dSwdU(NUM_VAR+2,NUM_VAR+2) += -kb[1];
//     dSwdU(NUM_VAR+2,NUM_VAR+3) += 0.5*M[2]*kf[1]*rho*c[4]/M[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
//     dSwdU(NUM_VAR+2,NUM_VAR+4) += M[2]*kf[1]/M[4]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);      
    
//     dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
//     dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];    
   
//     dSwdU(NUM_VAR+4,NUM_VAR) += 0.2*M[4]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*rho*c[1]/M[1]*pow(rho*c[1]/M[1],0.3)/M[0];
//     dSwdU(NUM_VAR+4,NUM_VAR+1) += rho*(0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)*c[1]*pow(rho*c[1]/M[1],0.5E-1)
// 		       *M[4]-0.25*kf[1]*c[4]*sqrt(rho*c[3]/M[3])*M[1])/(M[1]*M[1])/pow(rho*c_denom[1]/M[1],0.75);
//     dSwdU(NUM_VAR+4,NUM_VAR+2) += M[4]*kb[1]/M[2];
//     dSwdU(NUM_VAR+4,NUM_VAR+3) += -0.5*kf[1]*rho*c[4]/sqrt(rho*c_denom[3]/M[3])*pow(rho*c[1]/M[1],0.25)/M[3];
//     dSwdU(NUM_VAR+4,NUM_VAR+4) += -1.0*kf[1]*sqrt(rho*c[3]/M[3])*pow(rho*c[1]/M[1],0.25);

//     //this is a work around for the delta t calculation using an unesseccarily small value
//     //when CH4 & O2 -> ZERO
//     if( c_denom[0] != c[0] && CFL_flag ){
//       for(int i=0; i<num_react_species; i++){
// 	dSwdU(NUM_VAR+i,NUM_VAR+i)=ZERO; 
//       }
//     }
      
//    //take the minum value between laminar and turbulent concentrations
// //     if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
// // 	Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){
// //  }

//     break;
//     dSwdU(NUM_VAR,0) += M[0]*dkf_dT[0]*T/rho*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)+M[0]*
//       dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*
//       (-u*u-v*v+2.0*htot-2.0*Cp*p/rho/Rtot-2.0*c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)
//        -2.0*c[1]*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[3]*(h[3]-Cp*T*Rs[3]/Rtot-h[5]
//       +Cp*T*Rs[5]/Rtot)-2.0*c[4]*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot))/(Cp/Rtot-1.0)/2.0;

//     double top = pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1);

//     double sum_cs = c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot) + c[1]*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)
//       + c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot) + c[3]*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)
//       + c[4]*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot);

//     dSwdU(NUM_VAR,0) +=  (M[0]*dkf_dT[0]*T/rho)*top - (M[0]*dkf_dT[0]*top)*( (u*u+v*v)/2.0 - htot + Cp*T + sum_cs)/(rho*(Cp-Rtot));

//     dSwdU(NUM_VAR,1) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*u/(Cp/Rtot-1.0);
//     dSwdU(NUM_VAR,2) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*v/(Cp/Rtot-1.0);
//     dSwdU(NUM_VAR,3) += -M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)/(Cp/Rtot-1.0);

//     dSwdU(NUM_VAR,NUM_VAR) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)* pow(rho*c[1]/M[1],0.13E1)*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)
//       -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1);

//     dSwdU(NUM_VAR,NUM_VAR+1) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)
//       -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];

//     dSwdU(NUM_VAR,NUM_VAR+2) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);
//     dSwdU(NUM_VAR,NUM_VAR+3) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);
//     dSwdU(NUM_VAR,NUM_VAR+4) += M[0]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.13E1)*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);

//     MapleGenVar2 = M[1]*(0.15E1*dkf_dT[0]*T/rho*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.3*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],
// 0.13E1)*c[0]/M[0]-0.195E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],
// 0.3)*c[1]/M[1]+0.5*dkf_dT[1]*T*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M
// [1],0.25)-0.5*kf[1]*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)-0.25*kf[1]*rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)*c[3]/M[3]-0.125*kf[1]*rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho*c_denom
// [1]/M[1],0.75)*c[1]/M[1]-0.5*dkb_dT[1]*T*c[2]/M[2]+0.5*kb[1]*c[2]/M[2]);

//       MapleGenVar3 = -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(-u*u-v*v+2.0*htot-2.0*Cp*p/rho/Rtot-2.0*c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c
// [1]*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)
// -2.0*c[3]*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[4]*(h[4]-Cp*T*Rs[4]/Rtot-
// h[5]+Cp*T*Rs[5]/Rtot))/(Cp/Rtot-1.0)/2.0+0.3*M[1]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0]*c[0];

//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;
  
//       dSwdU(NUM_VAR+1,0) += MapleGenVar1-M[1]*(-0.195E1*kf[0]*pow(rho*c[0]/M[0],
// 0.2)*pow(rho*c[1]/M[1],0.3)*rho/M[1]-0.125*kf[1]*rho*rho*c[4]/M[4]*pow(rho
// *c[3]/M[3],0.5)/pow(rho*c_denom[1]/M[1],0.75)/M[1])*c[1]/rho-0.5*M[1]*kb[1]/M[2]
// *c[2]+0.25*M[1]*kf[1]*rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M
// [1],0.25)/M[3]*c[3]+0.5*M[1]*kf[1]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c
// [1]/M[1],0.25)*c[4]; 
      
    
//       dSwdU(NUM_VAR+1,1) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*
// 				   pow(rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*
// 				   pow(rho*c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*u/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+1,2) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*
// 				   pow(rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*
// 				   pow(rho*c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*v/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+1,3) += M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*
// 				  pow(rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*
// 				  pow(rho*c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])/(Cp/Rtot-1.0);
    
//       dSwdU(NUM_VAR+1,NUM_VAR) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/
// (Cp/Rtot-1.0)-0.3*M[1]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)
// /M[0];
//       dSwdU(NUM_VAR+1,NUM_VAR+1) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/
// (Cp/Rtot-1.0)+M[1]*(-0.195E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],
// 0.3)*rho/M[1]-0.125*kf[1]*rho*rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho
// *c_denom[1]/M[1],0.75)/M[1])/rho;
//       dSwdU(NUM_VAR+1,NUM_VAR+2) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/
// (Cp/Rtot-1.0)+0.5*M[1]*kb[1]/M[2];
//       dSwdU(NUM_VAR+1,NUM_VAR+3) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(
// rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/
// (Cp/Rtot-1.0)-0.25*M[1]*kf[1]*rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c
// [1]/M[1],0.25)/M[3];

//       dSwdU(NUM_VAR+1,NUM_VAR+4) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*
// 					   pow(rho*c[1]/M[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*
// 					   pow(rho*c[1]/M[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2])*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0) 
// 	                             -0.5*M[1]*kf[1]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25);
  

//     MapleGenVar2 = M[2]*(-dkf_dT[1]*T*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(
// rho*c[1]/M[1],0.25)+kf[1]*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M
// [1],0.25)+0.5*kf[1]*rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1]
// ,0.25)*c[3]/M[3]+0.25*kf[1]*rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho*c_denom
// [1]/M[1],0.75)*c[1]/M[1]+dkb_dT[1]*T*c[2]/M[2]-kb[1]*c[2]/M[2]);

//       MapleGenVar3 = -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(
// rho*c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*(-u*u-v*v+2.0*htot-2.0*Cp*p/rho/Rtot-2.0*c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c
// [1]*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)
// -2.0*c[3]*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[4]*(h[4]-Cp*T*Rs[4]/Rtot-
// h[5]+Cp*T*Rs[5]/Rtot))/(Cp/Rtot-1.0)/2.0-0.25*M[2]*kf[1]*rho*c[4]/M[4]*pow(rho*c
// [3]/M[3],0.5)/pow(rho*c_denom[1]/M[1],0.75)/M[1]*c[1];

//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;

//       dSwdU(NUM_VAR+2,0) += MapleGenVar1+kb[1]*c[2]-0.5*M[2]*kf[1]*rho*c[4]/M
// [4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)/M[3]*c[3]-M[2]*kf[1]/M
// [4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)*c[4];
//       dSwdU(NUM_VAR+2,1) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*u/(Cp/Rtot-1.0);

//       dSwdU(NUM_VAR+2,2) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*v/(Cp/Rtot-1.0);

//       dSwdU(NUM_VAR+2,3) += M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*
// c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+2,NUM_VAR) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp
// /Rtot-1.0);
//       dSwdU(NUM_VAR+2,NUM_VAR+1) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp
// /Rtot-1.0)+0.25*M[2]*kf[1]*rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho*c_denom[1]/
// M[1],0.75)/M[1];
     
//       dSwdU(NUM_VAR+2,NUM_VAR+2) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*
// 					   pow(rho*c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*
// 	(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)-kb[1];

//       dSwdU(NUM_VAR+2,NUM_VAR+3) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp
// /Rtot-1.0)+0.5*M[2]*kf[1]*rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M
// [1],0.25)/M[3];
//       dSwdU(NUM_VAR+2,NUM_VAR+4) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho
// *c[1]/M[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2])*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp
// /Rtot-1.0)+M[2]*kf[1]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25);

//       dSwdU(NUM_VAR+3,0) += -0.2E1*M[3]*dkf_dT[0]*T/rho*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)-0.1E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c
// [1]/M[1],0.13E1)*(-u*u-v*v+2.0*htot-2.0*Cp*p/rho/Rtot
// -2.0*c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[1]*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+
// Cp*T*Rs[5]/Rtot)-2.0*c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[3]*(h[3]-Cp*T*Rs
// [3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[4]*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot))/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,1) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*u/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,2) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*v/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,3) += 0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c
// [1]/M[1],0.13E1)/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,NUM_VAR) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)+0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
//       dSwdU(NUM_VAR+3,NUM_VAR+1) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)+0.26E1*M[3]*
// kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
//       dSwdU(NUM_VAR+3,NUM_VAR+2) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,NUM_VAR+3) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+3,NUM_VAR+4) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*
// c[1]/M[1],0.13E1)*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0);

//       MapleGenVar2 = M[4]*(-dkf_dT[0]*T/rho*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]
// /M[1],0.13E1)+0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)*c
// [0]/M[0]+0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)*c[1]/M
// [1]+dkf_dT[1]*T*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)-kf[1]*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)-0.5*kf[1]*
// rho*c[4]/M[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)*c[3]/M[3]-0.25*kf[1]
// *rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho*c_denom[1]/M[1],0.75)*c[1]/M[1]-
// dkb_dT[1]*T*c[2]/M[2]+kb[1]*c[2]/M[2]);

//       MapleGenVar3 = -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]
// /M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(-u*u-v*v+2.0*htot
// -2.0*Cp*p/rho/Rtot-2.0*c[0]*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[1]*(h[1]-Cp*
// T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[2]*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c
// [3]*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)-2.0*c[4]*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs
// [5]/Rtot))/(Cp/Rtot-1.0)/2.0-0.2*M[4]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/
// M[1],0.13E1)/M[0]*c[0];
//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;

//       dSwdU(NUM_VAR+4,0) += MapleGenVar1-M[4]*(0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)
// *pow(rho*c[1]/M[1],0.3)*rho/M[1]-0.25*kf[1]*rho*rho*c[4]/M[4]*pow(rho*c[3]
// /M[3],0.5)/pow(rho*c_denom[1]/M[1],0.75)/M[1])*c[1]/rho-M[4]*kb[1]/M[2]*c[2]+0.5
// *kf[1]*rho*c[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)/M[3]*c[3]+
// kf[1]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)*c[4];

//       dSwdU(NUM_VAR+4,1) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*u/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+4,2) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*v/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+4,3) += M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M
// [1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25
// )+dkb_dT[1]/Rtot*c[2]/M[2])/(Cp/Rtot-1.0);
//       dSwdU(NUM_VAR+4,NUM_VAR) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(h[0]-Cp*T*Rs[0]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)+0.2
// *M[4]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
   
//       dSwdU(NUM_VAR+4,NUM_VAR+1) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(h[1]-Cp*T*Rs[1]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)+M
// [4]*(0.13E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)*rho/M[1]
// -0.25*kf[1]*rho*rho*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)/pow(rho*c_denom[1]/M[1],
// 0.75)/M[1])/rho;
//       dSwdU(NUM_VAR+4,NUM_VAR+2) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(h[2]-Cp*T*Rs[2]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)+M
// [4]*kb[1]/M[2];
//       dSwdU(NUM_VAR+4,NUM_VAR+3) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(h[3]-Cp*T*Rs[3]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)-0.5
// *kf[1]*rho*c[4]/pow(rho*c_denom[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25)/M[3];
      
//       dSwdU(NUM_VAR+4,NUM_VAR+4) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/
// M[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]*pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],
// 0.25)+dkb_dT[1]/Rtot*c[2]/M[2])*(h[4]-Cp*T*Rs[4]/Rtot-h[5]+Cp*T*Rs[5]/Rtot)/(Cp/Rtot-1.0)-kf[1]
// *pow(rho*c[3]/M[3],0.5)*pow(rho*c[1]/M[1],0.25);


//     /******************** MODIFIED ****************************************/   
//     //with kf(T), kb(T)

//     kf[0] = reactions[0].kf(Temp);
//     kb[0] = ZERO;
//     kf[1] = reactions[1].kf(Temp);  
//     kb[1] = reactions[1].kb(Temp);
    
//     double dkf_dT[2], dkb_dT[2];

//     dkf_dT[0] = reactions[0].dkf_dT(Temp);
//     dkb_dT[0] = ZERO;
//     dkf_dT[1] = reactions[1].dkf_dT(Temp);
//     dkb_dT[1] = reactions[1].dkb_dT(Temp);

//     rho = W.rho;
//     double T = Temp;
//     double Cp = W.Cp();      
//     double u = W.v.x; 
//     double v = W.v.y; 
//     double p = W.p;      
//     double Rtot = W.Rtot();   
//     double htot = W.h();      
//     double h[6],R[6],CC[6],CC_denom[6];
//     double CP_Rtotm1 = Cp/T-ONE;
//     double MapleGenVar1, MapleGenVar2, MapleGenVar3;

//     double UNITS = 1e6;

//     for(int i=0; i<num_species; i++){    
//       M[i] = W.specdata[i].Mol_mass();
//       h[i] = (W.specdata[i].Enthalpy(Temp)+W.specdata[i].Heatofform());  
//       R[i] = (W.specdata[i].Rs()); 
//       CC[i] = rho*c[i]/M[i]/UNITS;              //W.SpecCon(i)/UNITS;      
//       CC_denom[i] = rho*c_denom[i]/M[i]/UNITS;
//     }    

// MOD WITH UNITS, CC, CC_denom
//       dSwdU(NUM_VAR,0) += M[0]*dkf_dT[0]*T/rho*pow(CC[0],0.2)*pow(CC[1],0.13E1)-M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(u*u+v*v-2.0*htot+2.0*Cp*p/rho/Rtot+2.0*c[0]*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[1]*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[2]*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[3]*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[4]*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot))/CP_Rtotm1/2.0;

//       dSwdU(NUM_VAR,1) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*u/CP_Rtotm1;
//       dSwdU(NUM_VAR,2) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*v/CP_Rtotm1;
//       dSwdU(NUM_VAR,3) += -M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)/CP_Rtotm1;
//       dSwdU(NUM_VAR,NUM_VAR) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.2*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/UNITS;

//       dSwdU(NUM_VAR,NUM_VAR+1) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.13E1*M[0]*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)/M[1]/UNITS;

//       dSwdU(NUM_VAR,NUM_VAR+2) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;     
//       dSwdU(NUM_VAR,NUM_VAR+3) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;
//       dSwdU(NUM_VAR,NUM_VAR+4) += M[0]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;

//       MapleGenVar2 = M[1]*(0.15E1*dkf_dT[0]*T/rho*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.3*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)*c[0]/M[0]/UNITS-0.195E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*c[1]/M[1]/UNITS+0.5*dkf_dT[1]*T*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-0.5*kf[1]*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-0.25*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)*c[3]/M[3]-0.125*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)*c[1]/M[1]-0.5*dkb_dT[1]*T*c[2]/M[2]/UNITS+0.5*kb[1]*c[2]/M[2]/UNITS);     

//       MapleGenVar3 = M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(u*u+v*v-2.0*htot+2.0*Cp*p/rho/Rtot+2.0*c[0]*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[1]*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[2]*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[3]*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[4]*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot))/CP_Rtotm1/2.0+0.3*M[1]*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/M[0]/UNITS*c[0];

//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;

//       dSwdU(NUM_VAR+1,0) += MapleGenVar1-M[1]*(-0.195E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*rho/M[1]/UNITS-0.125*kf[1]*rho*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1])*c[1]/rho-0.5*M[1]*kb[1]/M[2]/UNITS*c[2]+0.25*M[1]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3]*c[3]+0.5*M[1]*kf[1]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)*c[4];

//       dSwdU(NUM_VAR+1,1) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*u/CP_Rtotm1;

//       dSwdU(NUM_VAR+1,2) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*v/CP_Rtotm1;

//       dSwdU(NUM_VAR+1,3) += M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)/CP_Rtotm1;

//       dSwdU(NUM_VAR+1,NUM_VAR) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.3*M[1]*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/M[0]/UNITS;

//       dSwdU(NUM_VAR+1,NUM_VAR+1) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+M[1]*(-0.195E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*rho/M[1]/UNITS-0.125*kf[1]*rho*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1])/rho;

//       dSwdU(NUM_VAR+1,NUM_VAR+2) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.5*M[1]*kb[1]/M[2]/UNITS;

//       dSwdU(NUM_VAR+1,NUM_VAR+3) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.25*M[1]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3];

//       dSwdU(NUM_VAR+1, NUM_VAR+4) += -M[1]*(-0.15E1*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-0.5*dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.5*M[1]*kf[1]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25);

//       MapleGenVar2 = M[2]*(-dkf_dT[1]*T*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+kf[1]*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+0.5*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)*c[3]/M[3]+0.25*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)*c[1]/M[1]+dkb_dT[1]*T*c[2]/M[2]/UNITS-kb[1]*c[2]/M[2]/UNITS);

//       MapleGenVar3 = M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(u*u+v*v-2.0*htot+2.0*Cp*p/rho/Rtot+2.0*c[0]*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[1]*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[2]*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[3]*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[4]*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot))/CP_Rtotm1/2.0-0.25*M[2]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1]*c[1];

//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;

//       dSwdU(NUM_VAR+2,0) += MapleGenVar1+kb[1]/UNITS*c[2]-0.5*M[2]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3]*c[3]-M[2]*kf[1]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)*c[4];

//       dSwdU(NUM_VAR+2,1) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*u/CP_Rtotm1;
//       dSwdU(NUM_VAR+2,2) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*v/CP_Rtotm1;
//       dSwdU(NUM_VAR+2,3) += M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)/CP_Rtotm1;

//       dSwdU(NUM_VAR+2,NUM_VAR) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;

//       dSwdU(NUM_VAR+2,NUM_VAR+1) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.25*M[2]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1];

//       dSwdU(NUM_VAR+2,NUM_VAR+2) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-kb[1]/UNITS;

//       dSwdU(NUM_VAR+2,NUM_VAR+3) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.5*M[2]*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3];

//       dSwdU(NUM_VAR+2,NUM_VAR+4) += -M[2]*(dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+M[2]*kf[1]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25);

//       dSwdU(NUM_VAR+3,0) += -0.2E1*M[3]*dkf_dT[0]*T/rho*pow(CC[0],0.2)*pow(CC[1],0.13E1
// )+0.1E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(u*u+v*v-2.0*htot+2.0*Cp*p/rho/Rtot+2.0*c[0]*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[1]*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[2]*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[3]*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[4]*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot))/CP_Rtotm1;

//       dSwdU(NUM_VAR+3,1) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*u/CP_Rtotm1;
//       dSwdU(NUM_VAR+3,2) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*v/CP_Rtotm1;
//       dSwdU(NUM_VAR+3,3) += 0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)/CP_Rtotm1;

//       dSwdU(NUM_VAR+3,NUM_VAR) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.4*M[3]*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/M[0]/UNITS;

//       dSwdU(NUM_VAR+3,NUM_VAR+1) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.26E1*M[3]*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)/M[1]/UNITS;

//       dSwdU(NUM_VAR+3,NUM_VAR+2) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;
//       dSwdU(NUM_VAR+3,NUM_VAR+3) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;
//       dSwdU(NUM_VAR+3,NUM_VAR+4) += -0.2E1*M[3]*dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1;

//       MapleGenVar2 = M[4]*(-dkf_dT[0]*T/rho*pow(CC[0],0.2)*pow(CC[1],0.13E1)+0.2*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)*c[0]/M[0]/UNITS+0.13E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*c[1]/M[1]/UNITS+dkf_dT[1]*T*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-kf[1]*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)-0.5*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)*c[3]/M[3]-0.25*kf[1]*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)*c[1]/M[1]-dkb_dT[1]*T*c[2]/M[2]/UNITS+kb[1]*c[2]/M[2]/UNITS);

//       MapleGenVar3 = M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(u*u+v*v-2.0*htot+2.0*Cp*p/rho/Rtot+2.0*c[0]*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[1]*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[2]*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[3]*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)+2.0*c[4]*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot))/CP_Rtotm1/2.0-0.2*M[4]*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/M[0]/UNITS*c[0];

//       MapleGenVar1 = MapleGenVar2+MapleGenVar3;
//       dSwdU(NUM_VAR+4,0) += MapleGenVar1-M[4]*(0.13E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*rho/M[1]/UNITS-0.25*kf[1]*rho*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1])*c[1]/rho-M[4]*kb[1]/M[2]/UNITS*c[2]+0.5*kf[1]*rho*c[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3]*c[3]+kf[1]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)*c[4];

//       dSwdU(NUM_VAR+4,1) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*u/CP_Rtotm1;

//       dSwdU(NUM_VAR+4,2) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*v/CP_Rtotm1;

//       dSwdU(NUM_VAR+4,3) += M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)/CP_Rtotm1;

//       dSwdU(NUM_VAR+4,NUM_VAR) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[0]-Cp*T*R[0]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+0.2*M[4]*kf[0]/pow(CC_denom[0],0.8)*pow(CC[1],0.13E1)/M[0]/UNITS;

//       dSwdU(NUM_VAR+4,NUM_VAR+1) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[1]-Cp*T*R[1]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+M[4]*(0.13E1*kf[0]*pow(CC[0],0.2)*pow(CC[1],0.3)*rho/M[1]/UNITS-0.25*kf[1]*rho*rho*c[4]/M[4]/(UNITS*UNITS)*pow(CC[3],0.5)/pow(CC_denom[1],0.75)/M[1])/rho;

//       dSwdU(NUM_VAR+4,NUM_VAR+2) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[2]-Cp*T*R[2]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1+M[4]*kb[1]/M[2]/UNITS;

//       dSwdU(NUM_VAR+4,NUM_VAR+3) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[3]-Cp*T*R[3]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-0.5*kf[1]*rho*c[4]/(UNITS*UNITS)/pow(CC_denom[3],0.5)*pow(CC[1],0.25)/M[3];

//       dSwdU(NUM_VAR+4,NUM_VAR+4) += -M[4]*(dkf_dT[0]/rho/Rtot*pow(CC[0],0.2)*pow(CC[1],0.13E1)-dkf_dT[1]/Rtot*c[4]/M[4]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25)+dkb_dT[1]/Rtot*c[2]/M[2]/UNITS)*(h[4]-Cp*T*R[4]/Rtot-h[5]+Cp*T*R[5]/Rtot)/CP_Rtotm1-kf[1]/UNITS*pow(CC[3],0.5)*pow(CC[1],0.25);


/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables by using the complex step method
  
   dSwdU:  Matrix of source terms 
************************************************************************/
template<class SOLN_pSTATE>
inline void Reaction_set::Complex_Step_dSwdU(DenseMatrix &dSwdU, 
					     const SOLN_pSTATE &W) const {
  static const double TOL = 1.0E-15;
  static const double eps = 1.0E-100;
  static const cplx c_eps(0.0, eps);
  int NUM_VAR = W.NumVarSansSpecies();
  cplx PressDyne(TEN*W.p, ZERO); // pressure in dyne/cm^2
  cplx Temp(W.T(), ZERO);        // temperature in K 

  //assert(num_species == W.ns);

  // intitialize
  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND; //kg/mol -> g/mol 
    cc[i] = W.spec[i].c;                       //( W.spec[i].c < TOL )  ?  TOL : W.spec[i].c;
    //Wdot_perturbed[i] = 0.0;
  }

  //
  // compute the perturbed reaction rates moles/(cm**3*sec)
  //
  switch(reactset_flag) {
    
  //-------------------------------------------//
  //  15step CH4 mechanism based on GRI 2.11   //
  //-------------------------------------------//
  case CH4_15STEP_ARM2 :

    for(int j=0; j<num_react_species-1; j++){
      if ( W.spec[j].c > TOL ) {

	// perturb
	cc[j] += c_eps;
	
	// compute the perturbed reaction rates moles/(cm**3*sec)
	cplx15step211_(PressDyne, Temp, cc, rc);

	// set derivative
	for (int i=0; i < num_react_species-1; ++i) {
	  // in kg/m^3*s  ->  (g/mol)*(mol/cm^3*s)*1e3
	  rc[i] *= M[i]*THOUSAND;
	  dSwdU(NUM_VAR+i, NUM_VAR+j) = imag( rc[i]/eps )/W.rho;
	} // endfor


	// unperturbed value
	cc[j] -= c_eps;

      } // endif
    } // endfor 
    break;


  //-------------------------------------------//
  //  15step CH4 mechanism based on GRI 3
  //-------------------------------------------//
  case CH4_15STEP_ARM3 :

    for(int j=0; j<num_react_species-1; j++){
      if ( W.spec[j].c > TOL ) {
	
	// perturb
	cc[j] += c_eps;
	
	// compute the perturbed reaction rates moles/(cm**3*sec)
	cplx15step30_(PressDyne, Temp, cc, rc);

	// set derivative
	for (int i=0; i < num_react_species-1; ++i) {
	  // in kg/m^3*s  ->  (g/mol)*(mol/cm^3*s)*1e3
	  rc[i] *= M[i]*THOUSAND;
	  dSwdU(NUM_VAR+i, NUM_VAR+j) = imag( rc[i]/eps )/W.rho;
	} // endfor
	
	
	// unperturbed value
	cc[j] -= c_eps;

      } // endif
    } // endfor 
    break;
    
      
  default :
    cerr<<"\nReaction set "<<reactset_flag<<" is not valid"<<endl;
    exit(1);

  //-------------------------------------------//
  } // endswitch 
  //-------------------------------------------//
  

}  // end Complex_Step_dSwdU

/************************************************************************
  Calculates the Jacobian of the Chemical Source terms with respect
  to the conserved variables by using finite difference

  NOTE:  THis is a general way of doing it, i.e. calling omega(),
         but this is NOT EFFICIENT!!!!
  
   dSwdU:  Matrix of source terms 
************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE>
inline void Reaction_set::Finite_Difference_dSwdU(DenseMatrix &dSwdU, 
					     const SOLN_pSTATE &W,
					     const int Flow_Type) const {
  
  static const double b = numeric_limits<double>::epsilon();  // 1.0E-12;  
  int NUM_VAR = W.NumVarSansSpecies();
  SOLN_cSTATE U(W), S, S1, EPS; 

  // if this is one of those cases where num_species == num_reactions,
  // subtract one, cause dSwdU is (N-1)x(N-1)
  int ns=num_react_species;
  if (num_species == num_react_species) ns--;

  // get unperturbed reaction rate
  //S.Vacuum();
  omega<SOLN_pSTATE,SOLN_cSTATE>( S, W, Flow_Type );

  // Perturbation parameter (EPSILON)
  for (int j=0; j < ns; ++j) {
    EPS.rhospec[j].c = sqrt( b*fabs(U.rhospec[j].c) + b )/1000.0;
  }

  // initialize
  //S1 = S;

  //
  // compute finite differences
  //
  for (int j=0; j < ns; ++j) {
    if ( r[j] != 0.0 ) {

      // perturb
      // Because we are building the N-1 jacobian, when one species is 
      // perturbed, species N must be perturbed in the opposite direction
      // so that sum(mass_fracs)==1
      U.rhospec[j].c += EPS.rhospec[j].c;
      U.rhospec[num_species-1].c -= EPS.rhospec[j].c;

      // compute perturbed reaction rate
      omega<SOLN_pSTATE,SOLN_cSTATE>( S1, U.W(), Flow_Type );

      // set derivative
      for (int i=0; i < ns; ++i) {
	dSwdU(NUM_VAR+i, NUM_VAR+j) = TwoPointFiniteDifference(S.rhospec[i].c, 
							       S1.rhospec[i].c, 
							       EPS.rhospec[j].c);
      }

      // un-perturb
      U.rhospec[j].c -= EPS.rhospec[j].c;
      U.rhospec[num_species-1].c += EPS.rhospec[j].c;

    } // endif
  } // end for

} // end  Finite_Difference_dSwdU

#endif // _REACTIONS_INCLUDED
