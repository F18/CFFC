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

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif // _MATRIX_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

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
  double kf(const double &Temp,double &H2, double &O2, double &N2)const; //for H2&O2
  double kb(const double &Temp) const;
  
  template<class SOLN_pSTATE>
     double keq(const SOLN_pSTATE &W, const double &Temp) const;
  //Determine change in Gibbs free energy
  template<class SOLN_pSTATE>
     double deltaG(const SOLN_pSTATE &W) const;

  string react_name()const{ return react;}

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const React_data &W);
  friend istream& operator >> (istream &in_file,  React_data &W);
};

/*************** forward reaction coef (kf) ****************************/
inline double React_data::kf(const double &Temp) const{
  return A*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/*************** forward reaction kf for 2Step H2 ***********************/
inline double React_data::kf(const double &Temp,double &H2, double &O2, double &N2) const{
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


  return Astar*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/*************** backward reaction coef (kb) ****************************/
inline double React_data::kb(const double &Temp) const{
  return Ab*pow(Temp,n)*exp(-E/(R_UNIVERSAL*Temp));
}

/***************** equilibrium coef Keq ********************************/
template<class SOLN_pSTATE>
double React_data::keq(const SOLN_pSTATE &W,
                       const double& Temp) const{
   // nu_coef is the stoichiometric coef. sum (Eqn 13.25, 13.26 Anderson)
   // Kp or Keq has units of Pressure Pa( N/m^2) so need to change to
   // cgs units i.e *1e6
   return pow(R_UNIVERSAL*Temp, nu_coef)*exp(-deltaG(W)/(R_UNIVERSAL*Temp))*1e6;
}

/**********************************************************************
   Determine the Gibbs free energy change 
   deltaG = Sum(products) - Sum(reactants)
   
   Units:  J/mol
***********************************************************************/
template<class SOLN_pSTATE>
double React_data::deltaG(const SOLN_pSTATE &W) const{

  //need a general and more efficient system here !!!! SN
  // deltaG = sum (( product Gs) - (react Gs))

  //2STEP
  if(react == "H2O2_2step_1"){
    return TWO*W.Gibbs(2) - ( W.Gibbs(0) + W.Gibbs(1));
  } else if (react == "H2O2_2step_2"){
    return TWO*W.Gibbs(3) - ( TWO*W.Gibbs(2) + W.Gibbs(0));
 
    //8STEP
  } else if (react == "H2O2_8step_1"){
    return W.Gibbs(5) + W.Gibbs(0) - ( W.Gibbs(2) + W.Gibbs(1));
  } else if (react == "H2O2_8step_2"){
    return W.Gibbs(5) + W.Gibbs(2) - ( W.Gibbs(0) + W.Gibbs(3));
  } else if (react == "H2O2_8step_3"){
    return W.Gibbs(4) + W.Gibbs(2) - ( W.Gibbs(5) + W.Gibbs(3));
  } else if (react == "H2O2_8step_4"){
    return W.Gibbs(4) + W.Gibbs(0) - ( W.Gibbs(5) + W.Gibbs(5));
  } else if (react == "H2O2_8step_5"){
    return W.Gibbs(3) + W.Gibbs(4) - ( TWO*W.Gibbs(2) + W.Gibbs(4));
  } else if (react == "H2O2_8step_6"){
    return W.Gibbs(4) + W.Gibbs(4) - ( W.Gibbs(2) + W.Gibbs(5) + W.Gibbs(4));
  } else if (react == "H2O2_8step_7"){
    return W.Gibbs(5) + W.Gibbs(4) - ( W.Gibbs(2) + W.Gibbs(0) + W.Gibbs(4));
  } else if (react == "H2O2_8step_8"){
    return W.Gibbs(1) + W.Gibbs(4) - ( TWO*W.Gibbs(0) + W.Gibbs(4));
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
#define NO_REACTIONS 0

//Hardcoded Reaction Systems
#define CH4_1STEP 1
#define CH4_2STEP 2
#define H2O2_2STEP 3
#define H2O2_8STEP 4
#define CH4H2_2STEP 5

// User defined flag
#define USER 100

class Reaction_set{

private:  
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
  num_react_species=0; reactions = NULL; species = NULL;}
                      
  /******** Constructors *******************/
  //for hardcoded reactions
  void set_reactions(string &);
  //for user defined
  void set_species(string *, int);
  void set_reactions(int &,string*,double*,double*,double*);

  //Operator Overloading 
  Reaction_set& operator =(const Reaction_set &W);

  /* Input-output operators. */
  friend ostream& operator << (ostream &out_file, const Reaction_set &W);
  //friend istream& operator >> (istream &in_file,  Reaction_set &W);

  // time rate change of the species concentration 
  template<class SOLN_pSTATE, class SOLN_cSTATE>
     void omega( SOLN_cSTATE &U, 
                 const SOLN_pSTATE &W) const;
  template<class SOLN_pSTATE>
     void dSwdU(DenseMatrix &dSwdU, const SOLN_pSTATE &W,
                const bool &CFL_flag) const;

  
  void Deallocate();
  
  //destructor
  ~Reaction_set(){Deallocate();};
 
};


/**************** Destructor *******************************************/
inline void Reaction_set::Deallocate(){
  //deallocate memory
  if(reactions != NULL){
    delete[] reactions;
    reactions = NULL;
  }
  if(species != NULL){
    delete[] species;
    species = NULL;
  }
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

/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using the general law of mass action.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations
              of Euler3D_ThermallyPerfect_pState. (c_i*rho/M_i)   mol/m^3

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
template<class SOLN_pSTATE, class SOLN_cSTATE> 
void Reaction_set::omega(SOLN_cSTATE &U, 
                         const SOLN_pSTATE &W) const{

   double Temp = W.T();  //K
   double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3 
   double a,b;
   
   double *kf = new double[num_reactions];
   double *kb = new double[num_reactions];
   double *M = new double[num_react_species];
   double *c = new double[num_react_species];
   for(int i=0; i<num_react_species; i++){
      M[i] = W.specdata[i].Mol_mass()*THOUSAND;  //kg/mol -> g/mol
      c[i] = W.spec[i].c;                        //unitless
   }
   double ans = ZERO;

   switch(reactset_flag){
   case NO_REACTIONS:
      cerr<<"\n You shouldn't get here, NO_REACTIONS in Reaction_set::omeag(..)";
      exit(1);
      break;
 
   //ONE STEP CH4
   case CH4_1STEP: 
      //laminar case     
      a=1.0;// 0.2; 
      b=1.0;// 1.3;
      kf[0] = reactions[0].kf(Temp)*pow((W.SpecCon(0))/MILLION,a)*pow((W.SpecCon(1)/MILLION),b);              for(int index =0; index<num_react_species; index++){
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

  //clean up memory
  delete[] kf; delete[] kb;
  delete[] M; delete[] c;

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
template<class SOLN_pSTATE>
void Reaction_set::dSwdU(DenseMatrix &dSwdU, 
                         const SOLN_pSTATE &W, 
			 const bool &CFL_flag) const{

  /***************** Local Variables *****************/
   double Temp = W.T();
   double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3
   
   double *kf = new double[num_reactions];
   double *kb = new double[num_reactions];
   double *M = new double[num_react_species];
   double *c = new double[num_react_species];
   double *c_denom = new double[num_react_species];  //for denominator mass fracs
   double rate,tempvalue;
   double flag = 1.0;
   
   double tau_l, tau_t, tau_c, s, cm1;
   tau_l = ZERO;
   tau_c = ZERO;
   tau_t = ZERO;
   s = ZERO;
   cm1 = 4.0;

    
  for(int i=0; i<num_react_species; i++){
    M[i] = W.specdata[i].Mol_mass()*THOUSAND; //kg/mol -> g/mol 
    c[i] = W.spec[i].c;
    
    //For handling ~= ZERO mass fractions that appear in the denominator of dSwdU
    //by setting a lower tolerance allowed, and if below that set to tolerance 
    if( c[i] < TOLER){
      c_denom[i] = TOLER;
    }
    else {
      c_denom[i] = c[i];
    }
  }

  int NUM_VAR =  SOLN_pSTATE::NUM_VAR_3D - SOLN_pSTATE::ns;
   /*******************************************
   *  Reaction Mechanism Jacobians           *
   *                                         * 
   *******************************************/
  switch(reactset_flag){
  case NO_REACTIONS:
    //this case shouldn't be called
    cerr<<" Ummm, this dSwdU shouldn't be called for with NO_REACTIONS, please check your code.";
    break;
    
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
    // This is far from an elegant solution, but its easily created using maple.
    // It could probably be simplified for faster computation, but it works 
    // for now.

  //ONE STEP CH4
  case CH4_1STEP:
     
     kf[0] = reactions[0].kf(Temp);
     rate = kf[0]*pow((W.SpecCon(0))/MILLION,0.2)*pow((W.SpecCon(1)/MILLION),1.3);
     
     dSwdU(NUM_VAR,NUM_VAR) += -0.2*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1);
     dSwdU(NUM_VAR+1,NUM_VAR) += -0.4*M[1]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
     dSwdU(NUM_VAR+2,NUM_VAR) += 0.2*M[2]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
     dSwdU(NUM_VAR+3,NUM_VAR) += 0.4*M[3]*kf[0]/pow(rho*c_denom[0]/M[0],0.8)*pow(rho*c[1]/M[1],0.13E1)/M[0];
     
     dSwdU(NUM_VAR,NUM_VAR+1) += -0.13E1*M[0]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
     dSwdU(NUM_VAR+1,NUM_VAR+1) += -0.26E1*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3);
     dSwdU(NUM_VAR+2,NUM_VAR+1) += 0.13E1*M[2]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
     dSwdU(NUM_VAR+3,NUM_VAR+1) += 0.26E1*M[3]*kf[0]*pow(rho*c[0]/M[0],0.2)*pow(rho*c[1]/M[1],0.3)/M[1];
     //this is a work around for the delta t calculation using an unesseccarily small value
     //when c[0] -> ZERO
     if(c_denom[0] != c[0] && CFL_flag){ dSwdU(NUM_VAR,NUM_VAR)=ZERO; }
     break;

  //TWO STEP CH4   
  case CH4_2STEP:  
    
    //still some issues with units ??? ie.  dSwdU(6,6) ??
    //which is also on the diagonal so messes with CFL???
    kf[0] = reactions[0].kf(Temp);      
    kb[0] = ZERO;
    kf[1] = reactions[1].kf(Temp);  
    kb[1] = reactions[1].kb(Temp);

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

  /**************** Memory cleanup *****************/
  delete[] kf; delete[] kb;
  delete[] M; delete[] c;
  delete[] c_denom;

    
} //end dSwdU


#endif //REACTIONS
