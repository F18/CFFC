/****************** LESPremixedReactions.h **************************************
  This class defines the Reaction mechanisms for the 
  2D LES Premixed class.

  
          - based on Chem2D -> Reactions.h

***********************************************************************/
#ifndef _LESPREMIXED2D_REACTIONS_INCLUDED
#define _LESPREMIXED2D_REACTIONS_INCLUDED
 
class React_data;
class Reaction_set;

// Required C++ libraries
#include <iostream>
#include <string>
#include <cassert>

using namespace std;

#include "../Math/Math.h"
#include "../Math/Matrix.h"
#include "../Math/Tensor2D.h"
#include "../Physics/GasConstants.h"
#include "LESPremixed2DState.h"

#ifndef _COMPLEXIFY_INCLUDED
#include "../Math/Complexify/Complexify.h"
#endif // _COMPLEXIFY_INCLUDED


// List of subroutines not written in C/C++ 
extern "C" {
  // Reaction rates for CH4-15step mechanism based on GRI2.11
  void ckwyp15step211_( double& P, double& T, double* Y, double* Wdot);
  
  // Complexified version of the subroutine above
  void cplx15step211_( cplx& P, cplx& T, cplx* Y, cplx* Wdot);
  
  // Reaction rates for CH4-15step mechanism based on GRI3
  void ckwyp15step30_( double& P, double& T, double* Y, double* Wdot);

  // Complexified version of the subroutine above
  void cplx15step30_( cplx& P, cplx& T, cplx* Y, cplx* Wdot);
}


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
  double kf(const double &Temp, const double &H2, const  double &O2, const double &N2) const; //for H2&O2
  double kb(const double &Temp) const;
  double dkb_dT(const double &Temp) const;
  double keq(const LESPremixed2D_pState &W, const double &Temp) const;

  //Determine change in Gibbs free energy
  double deltaG(const LESPremixed2D_pState &W) const;

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
inline double React_data::keq(const LESPremixed2D_pState &W, const double& Temp) const{
  // nu_coef is the stoichiometric coef. sum (Eqn 13.25, 13.26 Anderson)
  // Kp or Keq has units of Pressure Pa( N/m^2) so need to change to 
  // cgs units i.e *1e6
  return pow(R_UNIVERSAL*Temp, nu_coef)*exp(-deltaG(W)/(R_UNIVERSAL*Temp))*1e6;
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
#define NO_REACTIONS 0
//Hardcoded Reaction Systems
#define CH4_1STEP          1
#define CH4_2STEP          2
#define H2O2_1STEP         3
#define H2O2_2STEP         4
#define H2O2_8STEP         5
#define CH4_15STEP_ARM2    6    // CH4 mechanisms based on GRI 2.11
#define CH4_15STEP_ARM3    7    // CH4 mechanisms based on GRI 3

// User defined flag
#define USER 100

class Reaction_set{

private:  

  //Temp vectors used during omega & dSwdU
  double *kf; 
  double *kb; 
  double *M; 
  double *c; 
  double *c_denom; 

protected:
public: 
  int reactset_flag;       //Reaction Set Flag
  React_data *reactions;   //each reaction rate data
  int num_reactions;       //number of reactions
  int num_species;         //number of total species
  int num_react_species;   //number of reacting species
  string *species;         //species used in reactions
  string Reaction_system;  //Reaction system name

  Reaction_set(){ reactset_flag=0; num_reactions=0; num_species=0; 
  num_react_species=0; reactions = NULL; species = NULL; 
  kf=NULL; kb=NULL; M=NULL; c=NULL; c_denom=NULL;}
                      
  /******** Constructors *******************/
  //for hardcoded reactions
  void set_reactions(string &);
  //for user defined
  void set_species(string *, int);
  void set_reactions(int &,string*,double*,double*,double*);

  //setup storage after num_reactions & num_species set.
  void set_storage(void){
    kf = new double[num_reactions];        
    kb = new double[num_reactions];
    M  = new double[num_react_species];
    c  = new double[num_react_species];
    c_denom = new double[num_react_species]; 
  }

  //Operator Overloading 
  Reaction_set& operator =(const Reaction_set &W);

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const Reaction_set &W);
  //friend istream& operator >> (istream &in_file,  Reaction_set &W);

  // time rate change of the species concentration 
  void omega(LESPremixed2D_cState &U, const LESPremixed2D_pState &W, const int Flow_Type ) const;

  //Jacobian ( flag true for cfl calc)
  void dSwdU(DenseMatrix &dSwdU,const LESPremixed2D_pState &W, const bool &CFL_flag, const int Flow_Type, const int Solver_type) const;

  void Finite_Difference_dSwdU(DenseMatrix &dSwdU, const LESPremixed2D_pState &Wlocal, const int Flow_Type) const;
  
  void Complex_Step_dSwdU(DenseMatrix &dSwdU, const LESPremixed2D_pState &Wlocal) const;

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
}


/***************** Assignment ****************************************/
inline Reaction_set& Reaction_set::operator =(const Reaction_set &W){
  //self assignment protection
  if( this != &W){   
    string temp = W.Reaction_system;
    //copy assignment
    set_reactions(temp);
  }
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

#endif // _LESPREMIXED2D_REACTIONS_INCLUDED
