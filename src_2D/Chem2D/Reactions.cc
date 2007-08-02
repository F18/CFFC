/****************** Reactions.cc ************************************
  This class defines the Reaction class and for the 
  2D Axisymmertric Navier-Stokes with Multiple Species solution.

***********************************************************************/
#ifndef _REACTIONS_INCLUDED
#include "Reactions.h"
#endif // _REACTIONS_INCLUDED

/***********************************************************************
  Setup reaction data using hardcoded reaction mechanisms.
***********************************************************************/
void Reaction_set::set_reactions(string &react){

  Deallocate();

  Reaction_system = react;
  //no reactions
  if( react == "NO_REACTIONS"){
    reactset_flag = NO_REACTIONS;
    reactions = NULL;
    num_reactions = 0;
  }

  /****************** METHANE ***********************************/   
  // 1step CH4 mechanism
  else if( react == "CH4_1STEP" ){
    //flag for hadcoded
    reactset_flag =  CH4_1STEP;
    //set number of reactions & species
    num_reactions = 1;
    num_species = 5;
    num_react_species = 4;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients based on kov formulation
    reactions[0].set_data("CH4_1step",3.3*6.7e12,0.0,(48400.0*CAL_TO_JOULE),0);
        
    //setup species list 
    species = new string[num_species];
    species[0] = "CH4";
    species[1] = "O2";
    species[2] = "CO2";
    species[3] = "H2O";
    species[4] = "N2";
  }
   
  // 2step CH4 mechanism 
  else if(react == "CH4_2STEP"){
    //flag for hardcoded
    reactset_flag = CH4_2STEP;
    //set number of reactions & species
    num_reactions = 2;
    num_species = 6;
    num_react_species = 5;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients
    reactions[0].set_data("CH4_2step_1",2.6*1.443e13,0.0,(48400.0*CAL_TO_JOULE),0);
    reactions[1].set_data("CH4_2step_2",pow(10.0,14.6),5.0e8,0.0,(40000.0*CAL_TO_JOULE),0);
    
    //setup species list 
    species = new string[num_species];
    species[0] = "CH4";
    species[1] = "O2";
    species[2] = "CO2";
    species[3] = "H2O";
    species[4] = "CO"; 
    species[5] = "N2";
  }

  /****************** HYDROGEN ***********************************/   
  //2step H2 & O2
  else if(react == "H2O2_2STEP"){
    //flag for hardcoded
    reactset_flag = H2O2_2STEP;
    //set number of reactions & species
    num_reactions = 2;
    num_species = 5;
    num_react_species = 4;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients
    reactions[0].set_data("H2O2_2step_1",1.0,-10.0,4865.0*CAL_TO_JOULE,0);
    reactions[1].set_data("H2O2_2step_2",1.0,-13.0,42000.0*CAL_TO_JOULE,1);
  
    //setup species list 
    species = new string[num_species];
    species[0] = "H2";
    species[1] = "O2";
    species[2] = "OH";
    species[3] = "H2O"; 
    species[4] = "N2";
  }
  
  else if(react == "H2O2_8STEP"){
    //8Step H2 & O2
    //flag for hardcoded
    reactset_flag = H2O2_8STEP;
    //set number of reactions & species
    num_reactions = 8;
    num_species = 7;
    num_react_species = 6;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients
    reactions[0].set_data("H2O2_8step_1",2.2e14,0   ,16800.0*CAL_TO_JOULE,0);
    reactions[1].set_data("H2O2_8step_2",1.8e10,1.0 ,8900.0*CAL_TO_JOULE ,0);
    reactions[2].set_data("H2O2_8step_3",2.2e13,0   ,5150.0*CAL_TO_JOULE ,0);
    reactions[3].set_data("H2O2_8step_4",6.3e12,0   ,1090.0*CAL_TO_JOULE   ,0);
    reactions[4].set_data("H2O2_8step_5",6.6e17,-1.0,0.0*CAL_TO_JOULE    ,1);
    reactions[5].set_data("H2O2_8step_6",2.2e22,-2.0,0.0*CAL_TO_JOULE    ,1);
    reactions[6].set_data("H2O2_8step_7",6.0e16,-0.6,0.0*CAL_TO_JOULE    ,1);
    reactions[7].set_data("H2O2_8step_8",6.0e13,0   ,-1000.0*CAL_TO_JOULE  ,1);

    species = new string[num_species];
    species[0] = "O";
    species[1] = "O2";
    species[2] = "H";
    species[3] = "H2"; 
    species[4] = "H2O";
    species[5] = "OH";
    species[6] = "N2";
  }
  
  // else error  
  else{
    cerr<<"\nReaction set "<<react<<" is not valid"<<endl;
    exit(1);
  }

} //end of set_reactions


/***********************************************************************
  Setup user defined species.  Basically copy it to a perminant
  storage location so pointers can be used.
***********************************************************************/
void Reaction_set::set_species(string *spec, int num){
  Deallocate();
  num_species=num;
  species = new string[num];
  for(int i=0; i<num; i++){
    species[i] = spec[i];
  }
}

/***********************************************************************
  Setup user defined mechanisms.

  NOT FINISHED!!!!
***********************************************************************/
void Reaction_set::set_reactions(int &num_react,string* name, double* A,
			   double* b, double* E){  
  
//   //flag for user defined
//   reactset_flag=USER;
//   //number of reactions in set
//   num_reactions = num_react;   
//   //allocate memory 
//   reactions = new React_data[num_react];    
//   //get each reaction set data
//   for(int i=0; i<num_react; i++){
//     reactions[i].set_data(name[i],A[i],b[i],E[i]);
//   }
  cout<<"\n USER REACTIONS NOT AVAILABE YET ";
  cout<<"\n Please send money ";
  exit(1);

} //end of set_reactions

/**********************************************************************
   Determine the Gibbs free energy change 
   deltaG = Sum(products) - Sum(reactants)
   
   Units:  J/mol
***********************************************************************/
double React_data::deltaG(const Chem2D_pState &W) const{

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

/************************************************************************
  Calculates the concentration time rate of change of species from
  primitive state W using the general law of mass action.
  U is the conserved state container for passing back the 
  source terms. ie. U.rhospec[i].c 

  W.SpecCon:  is the  species mass fractions concentrations
              of Chem2D_pState. (c_i*rho/M_i)   mol/m^3

  Return units are  kg/m^3*s ie. rho*omega (kg/m^3)*(1/s)

************************************************************************/
void Reaction_set::omega(Chem2D_cState &U, const Chem2D_pState &W,  const int Flow_Type ) const{
 
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
    //---------------------------------//
    //------ Hardcoded ----------------//
    //---------------------------------//
  case CH4_1STEP: 
    //laminar case     
    a = 0.2; 
    b = 1.3;
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
    //turbulent case
    if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA ||
	Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON) {
   
      double rhoc_minum =  U.rhospec[0].c; 
      for(int index =1; index<num_react_species; index++){
	if (U.rhospec[index].c < rhoc_minum)
 	  {
	    rhoc_minum =U.rhospec[index].c;
	  }
      } // Find the deficit species in terms of concentration
      for(int index =0; index<num_react_species; index++){
	switch(index) {
	case 0 : //CH4
	  ans = -W.Coeff_edm*U.rhoomega*rhoc_minum/U.rho;
	  break;
	case 1 : //O2
	  ans = -HALF*W.Coeff_edm*U.rhoomega*rhoc_minum/U.rho;
	  break;
	case 2 : //CO2
	  ans = W.Coeff_edm*U.rhoomega*rhoc_minum/U.rho;
	  break;
	case 3 : //H2O
	  ans = TWO*W.Coeff_edm*U.rhoomega*rhoc_minum/U.rho;
	  break;
	};
	//Use the minmum reaction rate between laminar and turbulent 
	if(index<2){
	  U.rhospec[index].c = max(ans, U.rhospec[index].c);  
	}else{
	  U.rhospec[index].c = min(ans, U.rhospec[index].c); 
	  //This will be useful for chemistry source Jacobian because
	  //the limitting case (finite rate chemisty and turbulent mixing) is the one
	  //controlling the reaction process ... 
	  //For forward reaction scheme actually only the fuel and oxidizer concentration (max) 
	  //are useful, for products (max) which might be reconsidered when considering reversible 
	//reaction scheme 
	}	 
	
      }
      
    }//turbulent case  
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
void Reaction_set::dSwdU(DenseMatrix &dSwdU, const Chem2D_pState &W, 
			 const bool &CFL_flag, const int Flow_Type) const{

  /***************** Local Variables *****************/
  double Temp = W.T();
  double rho= W.rho/THOUSAND; //kg/m^3 -> g/cm^3
  
  double *kf = new double[num_reactions];
  double *kb = new double[num_reactions];
  double *M = new double[num_react_species];
  double *c = new double[num_react_species];
  double *c_denom = new double[num_react_species];  //for denominator mass fracs

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

  int NUM_VAR = NUM_CHEM2D_VAR_SANS_SPECIES;
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
    // It could probably be simplified for faster computation, but it works 
    // for now.

  case CH4_1STEP:
    kf[0] = reactions[0].kf(Temp); 

    // One forward reaction ... so only depending on the fuel and oxidizer
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
    
    //take the minum value between laminar and turbulent concentrations
    if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
 	Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){
      
      double c_min = 0.0;
      //fuel and oxidizer [0],[1]
      for(int i=0; i<2; i++){
	c_min = W.spec[0].c;
	if(W.spec[i].c<c_min){
	  c_min =W.spec[i].c;
	}
	
      }//find out the deficient species in fuel and oxidizer (edm)
      for(int i=0; i<2; i++){
	dSwdU(NUM_VAR+i,0) = -W.Coeff_edm *c_min*W.omega;
	dSwdU(NUM_VAR+i,5) = -W.Coeff_edm*W.rho*c_min;
      }
 
      //The limiting case has been determined in Sw(...)
      Chem2D_cState Sw_LAMINAR;
      Chem2D_cState Sw_TURBULENT;
      Sw_LAMINAR = W.Sw(W.React.reactset_flag,FLOWTYPE_LAMINAR);
      Sw_TURBULENT= W.Sw(W.React.reactset_flag,FLOWTYPE_TURBULENT_RANS_K_OMEGA);
      
      if (Sw_LAMINAR[0]<Sw_TURBULENT[0]){
	dSwdU(NUM_VAR,NUM_VAR) =  -W.rho* W.Coeff_edm*W.omega;
	dSwdU(NUM_VAR+1,NUM_VAR) = -TWO*dSwdU(NUM_VAR,NUM_VAR);
	dSwdU(NUM_VAR+2,NUM_VAR) = dSwdU(NUM_VAR,NUM_VAR);
	dSwdU(NUM_VAR+3,NUM_VAR) = TWO*dSwdU(NUM_VAR,NUM_VAR);
	dSwdU(NUM_VAR,NUM_VAR+1) =  -W.rho* W.Coeff_edm*W.omega;
	dSwdU(NUM_VAR+1,NUM_VAR+1) = -TWO*dSwdU(NUM_VAR,NUM_VAR+1) ;
	dSwdU(NUM_VAR+2,NUM_VAR+1) = dSwdU(NUM_VAR,NUM_VAR);
	dSwdU(NUM_VAR+3,NUM_VAR+1) = TWO*dSwdU(NUM_VAR,NUM_VAR);
	
      }
      
    }//end of turbulent case
    
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

    //this is a work around for the delta t calculation using an unesseccarily small value
    //when CH4 & O2 -> ZERO
    if( c_denom[0] != c[0] && CFL_flag ){
      for(int i=0; i<num_react_species; i++){
	dSwdU(NUM_VAR+i,NUM_VAR+i)=ZERO; 
      }
    }
      
   //take the minum value between laminar and turbulent concentrations
//     if (Flow_Type == FLOWTYPE_TURBULENT_RANS_K_EPSILON ||
// 	Flow_Type == FLOWTYPE_TURBULENT_RANS_K_OMEGA){

//  }

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
