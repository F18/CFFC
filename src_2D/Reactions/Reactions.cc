/****************** Reactions.cc ************************************
  This class defines the Reaction class and for the 
  2D Axisymmertric Navier-Stokes with Multiple Species solution.

***********************************************************************/
#ifndef _REACTIONS_INCLUDED
#include "Reactions.h"
#endif // _REACTIONS_INCLUDED

//
// Cantera Libraries
//
#ifdef _CANTERA_VERSION
#include <cantera/equilibrium.h>  // canteraequilibrium solver
#endif //_CANTERA_VERSION

/***********************************************************************
  Static object initialization.
***********************************************************************/
// indexes of specific species
int Reaction_set::iCO  = -1;
int Reaction_set::iCO2 = -1;
int Reaction_set::iH2O = -1;
int Reaction_set::iO2  = -1;

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

    // set static indexes (-1 if does not exist)
    iCO  = -1;
    iH2O = -1;
    iCO2 = -1;
    iO2  = -1;
  }

  /****************** METHANE ***********************************/   
  //-------------------------------------------//
  // 1step CH4 mechanism
  //-------------------------------------------//
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
    //reactions[0].set_data("CH4_1step",3.3*6.7e12,0.0,(48400.0*CAL_TO_JOULE),0);    for 0.2,1,3
    reactions[0].set_data("CH4_1step",4.5*2.4e16,0.0,(48400.0*CAL_TO_JOULE),0);     //for 1.0,1.0
      
    //setup species list 
    species = new string[num_species];
    species[0] = "CH4";
    species[1] = "O2";
    species[2] = "CO2";
    species[3] = "H2O";
    species[4] = "N2";

    // set static indexes (-1 if does not exist)
    iCO  = -1;
    iH2O = 3;
    iCO2 = 2;
    iO2  = 1;
  }
   
  //-------------------------------------------//
  // 2step CH4 mechanism 
  //-------------------------------------------//
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
 
    // set static indexes (-1 if does not exist)
    iCO  = 4;
    iH2O = 3;
    iCO2 = 2;
    iO2  = 1;
 }

  // 1step C3H8 mechanism
  else if( react == "C3H8_1STEP" ){
    //flag for hadcoded
    reactset_flag =  C3H8_1STEP;
    //set number of reactions & species
    num_reactions = 1;
    num_species = 5;
    num_react_species = 4;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients based on kov formulation
    //reactions[0].set_data("CH4_1step",3.3*6.7e12,0.0,(48400.0*CAL_TO_JOULE),0);    for 0.2,1,3
    reactions[0].set_data("C3H8_1step",4.5*2.4e16,0.0,(48400.0*CAL_TO_JOULE),0);     //for 1.0,1.0
      
    //setup species list 
    species = new string[num_species];
    species[0] = "C3H8";
    species[1] = "O2";
    species[2] = "CO2";
    species[3] = "H2O";
    species[4] = "N2";
  }

  //-------------------------------------------//
  //  15step CH4 mechanism based on GRI 2.11   //
  //-------------------------------------------//
  else if(react == "CH4_15STEP_ARM2") {
    //flag for hardcoded
    reactset_flag = CH4_15STEP_ARM2;
    //set number of reactions & species
    num_reactions = 0;
    num_species = 19;
    num_react_species = 19; // All needed by the fortran subroutine

    //setup species list 
    species = new string[num_species];
    species[0]  = "H2";               
    species[1]  = "H";               
    species[2]  = "O2";              
    species[3]  = "OH";               
    species[4]  = "H2O";              
    species[5]  = "HO2";              
    species[6]  = "H2O2";            
    species[7]  = "CH3";            
    species[8]  = "CH4";             
    species[9]  = "CO";              
    species[10] = "CO2";              
    species[11] = "CH2O";            
    species[12] = "C2H2";            
    species[13] = "C2H4";             
    species[14] = "C2H6";            
    species[15] = "NH3";              
    species[16] = "NO";               
    species[17] = "HCN";              
    species[18] = "N2";

    // set static indexes (-1 if does not exist)
    iCO  = 9;
    iH2O = 4;
    iCO2 = 10;
    iO2  = 2;

  /***********************
    Reactions:

    (1)    H +  0.50 O2 = OH 
    (2)    H2 +  0.50 O2 = H + OH                                                
    (3)    HO2 =  0.50 O2 + OH                                         
    (4)    0.50 O2 + H2O2 = OH + HO2
    (5)    0.50 O2 +  0.50 C2H2 = H + CO 
    (6)    CH3 + CO + C2H4 =  0.50 O2 + CH4 +  1.50 C2H2
    (7)    0.50 O2 +  2CH3 = H2 + CH4 + CO 
    (8)    0.50 O2 + CH3 = H + CH2O                                                
    (9)    0.50 O2 + CH4 = OH + CH3
    (10)   0.50 O2 + CO = CO2
    (11)   0.50 O2 + C2H6 = CH4 + CH2O
    (12)   H + OH = H2O 
    (13)   H + CH4 + NO + HCN =  0.50 O2 +  2CH3 + N2  
    (14)   H +  0.50 O2 + CH4 + HCN =  2CH3 + NO
    (15)   0.50 O2 + CH4 + NH3 + HCN = H2O +  2CH3 + N2
                   
  ***********************/
  }

  //-------------------------------------------//
  //  15step CH4 mechanism based on GRI 3
  //-------------------------------------------//
  else if(react == "CH4_15STEP_ARM3") {
    //flag for hardcoded
    reactset_flag = CH4_15STEP_ARM3;
    //set number of reactions & species
    num_reactions = 0;
    num_species = 19;
    num_react_species = 19;  // All needed by the fortran subroutine

    //setup species list 
    species = new string[num_species];
    species[0]  = "H2";               
    species[1]  = "H";               
    species[2]  = "O2";              
    species[3]  = "OH";               
    species[4]  = "H2O";              
    species[5]  = "HO2";              
    species[6]  = "H2O2";            
    species[7]  = "CH3";            
    species[8]  = "CH4";             
    species[9]  = "CO";              
    species[10] = "CO2";              
    species[11] = "CH2O";            
    species[12] = "C2H2";            
    species[13] = "C2H4";             
    species[14] = "C2H6";            
    species[15] = "NH3";              
    species[16] = "NO";               
    species[17] = "HCN";              
    species[18] = "N2";

    // set static indexes (-1 if does not exist)
    iCO  = 9;
    iH2O = 4;
    iCO2 = 10;
    iO2  = 2;

  /***********************
    Reactions:

    ( 1)   2H +  2OH =  2H2 + O2                                             
    ( 2)   2H = H2                                                           
    ( 3)   H + HO2 = H2 + O2                                                 
    ( 4)   H + H2O2 = H2 + HO2                                                
    ( 5)   OH + CH3 = H2 + CH2O                                               
    ( 6)   H + CH4 = H2 + CH3                                                 
    ( 7)   H + OH + CO = H2 + CO2                                             
    ( 8)   CH2O = H2 + CO                                                     
    ( 9)   O2 + C2H2 = H2 +  2CO                                              
    (10)   OH + C2H4 = H2 + CH3 + CO                                          
    (11)   C2H6 = H2 + C2H4                                                   
    (12)   H + OH = H2O                                                       
    (13)   2NO = O2 + N2                                                     
    (14)   H2 + CO + NO = H + O2 + HCN                                        
    (15)   3H + H2O + NH3 =  4H2 + NO                   
  ***********************/
  }


  /****************** HYDROGEN ***********************************/   
  //-------------------------------------------//
  // 1step H2 & O2
  //-------------------------------------------//
  else if(react == "H2O2_1STEP"){    //UNTESTED !!!!!!!!
    //flag for hardcoded
    reactset_flag = H2O2_1STEP;
    //set number of reactions & species
    num_reactions = 2;
    num_species = 4;
    num_react_species = 3;
    //allocate memory 
    reactions = new React_data[num_reactions];
    //set reaction coefficients
    reactions[0].set_data("H2O2_1step_1",4.0e5,0.0,7.0*R_UNIVERSAL*1000.0,0);
  
    //setup species list 
    species = new string[num_species];
    species[0] = "H2";
    species[1] = "O2";
    species[2] = "H2O"; 
    species[3] = "N2";

    // set static indexes (-1 if does not exist)
    iCO  = -1;
    iH2O = 2;
    iCO2 = -1;
    iO2  = 1;
 }

  //-------------------------------------------//
  // 2step H2 & O2
  //-------------------------------------------//
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

    // set static indexes (-1 if does not exist)
    iCO  = -1;
    iH2O = 3;
    iCO2 = -1;
    iO2  = 1;
  }
  
  //-------------------------------------------//
  // 8step H2 & O2
  //-------------------------------------------//
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

    // set static indexes (-1 if does not exist)
    iCO  = -1;
    iH2O = 4;
    iCO2 = -1;
    iO2  = 1;
  }
  
  // else error  
  else{
    cerr<<"\nReaction set "<<react<<" is not valid"<<endl;
    exit(1);
  }

  //SETUP TEMP STORAGE 
  if( react != "NO_REACTIONS"){
    set_storage();
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

/**************************************************************
  Get the index corrensponding to the species name.  
  If it is not found, return -1;                      
***************************************************************/
int Reaction_set::SpeciesIndex( const string& name ) const {
  int index(-1);
  for(int k=0; k<num_species; k++) {
    if (species[k]==name) { index = k; break; }
  }
  return index;
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


/***********************************************************************
  Load Cantera mechanisms.  A cantera IdealGasMix object is created
  from the user supplied mechanism file and specific mechanism name.
  Basically, Cantera has it's own parser for either chemkin (*.ck), 
  or cantera (*.cti) formated mechanism files.  The IdealGasMix
  stores the mechanism data (such as reaction constants, high/low 
  pressure limits) and has member functions that allow for the 
  computation of reaction rates.  
***********************************************************************/
void Reaction_set::ct_load_mechanism(string &mechanism_file_name, 
				     string &mechanism_name) 
{

// _CANTERA_VERSION flag set
#ifdef _CANTERA_VERSION
 
 // make sure all unused parameters are null
  Deallocate();

  //flag for cantera
  reactset_flag = CANTERA;
  
  //create a new ideal gas mixture class
  try {
    ct_gas = new IdealGasMix(mechanism_file_name, mechanism_name);
  }
  catch (CanteraError) {
    Cantera::showErrors();
  }

  //get the number of reactions and species
  num_species = ct_gas->nSpecies();
  num_react_species = num_species;
  num_reactions = ct_gas->nReactions();

  //set the species names 
  species = new string[num_species];
  for(int index =0; index<num_species; index++){
     species[index] = ct_gas->speciesName(index);
  }

  //set the reaction system name
  ct_mech_name = ct_gas->name();
  ct_mech_file = mechanism_file_name;
  Reaction_system = "CANTERA";

  // set static indexes (-1 returbed if does not exist)
  iCO  = ct_gas->speciesIndex("CO");
  iH2O = ct_gas->speciesIndex("H2O");
  iCO2 = ct_gas->speciesIndex("CO2");
  iO2  = ct_gas->speciesIndex("O2");

  // allocate some temporary storage
  set_storage();



// _CANTERA_VERSION flag not set
#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION

} //end of ct_load_mechanism

/***********************************************************************
  Use cantera to parse the input mass fraction string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have 0 mass fractions.  Cantera
  also normalizes the mass fractions to sum to unity.  Returns them
  in an array.
***********************************************************************/
void Reaction_set::ct_parse_mass_string(const string& massFracStr, 
					double* massFracs) {

#ifdef _CANTERA_VERSION

  //
  // for a string of 'space' delimited mass fractions
  //
  if ( massFracStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = massFracStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = massFracStr.find_first_of(delimiters, lastPos);

    // parse the string and set mass fractions
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      massFracs[index] = atof(massFracStr.substr(lastPos, pos - lastPos).c_str());
      massFracs[index] = max(massFracs[index], ZERO);
      lastPos = massFracStr.find_first_not_of(delimiters, pos);
      pos = massFracStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

    // normalize to unity
    double sum(0.0);
    for(int index=0; index<num_species; index++) sum += massFracs[index];
    if (sum>TOLER) for(int index=0; index<num_species; index++) massFracs[index] /= sum;

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // set mass fractions and make sure they sum to unity
    ct_gas->setMassFractionsByName(massFracStr);
    for(int index =0; index<num_species; index++)
      massFracs[index] = ct_gas->massFraction(index);

  } // endif

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";  
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_parse_mass_string

/***********************************************************************
  Use cantera to parse the input mole fraction string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have 0 mole fractions.  Cantera
  also normalizes the mole fractions to sum to unity.  Returns them in
  an array.
***********************************************************************/
void Reaction_set::ct_parse_mole_string(const string& moleFracStr, 
					double* moleFracs) {

#ifdef _CANTERA_VERSION

  //
  // for a string of 'space' delimited mole fractions
  //
  if ( moleFracStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = moleFracStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = moleFracStr.find_first_of(delimiters, lastPos);

    // parse the string and set molar fractions
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      moleFracs[index] = atof(moleFracStr.substr(lastPos, pos - lastPos).c_str());
      moleFracs[index] = max(moleFracs[index], ZERO);
      lastPos = moleFracStr.find_first_not_of(delimiters, pos);
      pos = moleFracStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

    // normalize to unity
    double sum(0.0);
    for(int index=0; index<num_species; index++) sum += moleFracs[index];
    if (sum>TOLER) for(int index=0; index<num_species; index++) moleFracs[index] /= sum;

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // set mole fractions and make sure they sum to unity
    ct_gas->setMoleFractionsByName(moleFracStr);
    for(int index =0; index<num_species; index++)
      moleFracs[index] = ct_gas->moleFraction(index);

  } // endif

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";  
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_parse_mass_string


/***********************************************************************
  Use cantera to parse the input schmidt number string of the form
      CH4:0.5, O2:0.5
  All other species will be assumed to have unity Schmidt number.  
  Returns them in an array.
***********************************************************************/
void Reaction_set::ct_parse_schmidt_string( const string& schmidtStr, 
					    double* schmidt) {

#ifdef _CANTERA_VERSION

  //
  // for a string of 'space' delimited schmidt numbers
  //
  if ( schmidtStr.find(":",0) == string::npos ) {
    
    // get first position
    string delimiters = " ";
    string::size_type lastPos = schmidtStr.find_first_not_of(delimiters, 0);
    string::size_type pos     = schmidtStr.find_first_of(delimiters, lastPos);

    // parse the string and set schmidt numbers
    int index = 0;
    while (string::npos != pos || string::npos != lastPos) {
      schmidt[index] = atof(schmidtStr.substr(lastPos, pos - lastPos).c_str());
      if(schmidt[index]<0.0) schmidt[index] = ONE;
      lastPos = schmidtStr.find_first_not_of(delimiters, pos);
      pos = schmidtStr.find_first_of(delimiters, lastPos);
      index++;
    } // endwhile

  //
  // for a composition map (eg. CH4:0.5, O2:0.5)
  //
  } else {

    // declares
    compositionMap xx;
    int kk = ct_gas->nSpecies();
    double s;
    
    // build composition map and initialize
    for (int k = 0; k < kk; k++) xx[ct_gas->speciesName(k)] = -1.0;
    
    // parse map
    parseCompString(schmidtStr, xx);
    
    // set schmidt numbers
    for (int k = 0; k < kk; k++) { 
      s = xx[ct_gas->speciesName(k)];
      if (s > ZERO) schmidt[k] = s;
      else schmidt[k] = ONE;
    }
    
  } // endif

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_parse_mass_string

/***********************************************************************
  Use cantera to return the species index.  Exits in error if an
  unidentified species is requested.
***********************************************************************/
int Reaction_set::ct_get_species_index(const string &sp) {

#ifdef _CANTERA_VERSION

  // get index
  int index = ct_gas->speciesIndex(sp);

  // if index==-1, species not found
  if (index<0) {
    cerr<<"\n Reaction_set::ct_get_species_index() - Index of unkown species, "
	<<sp<<", requested.\n";
    exit(-1);
  } // endif

  // return index
  return index;


#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_get_species_index


/***********************************************************************

  This function is required as a workarround due to conflicts between
  ::DenseMatrix and Cantera::DenseMatrix.  Instead of inserting '::'
  everywhere in front of DenseMatrix, I just move the header file into
  this include.  The problem only comes about when including 
  <cantera/equilibrium.h> which defines the function 
  'Cantera::equilibrate()'.

  All this function does is equilibrate the fresh gas mixture defined
  in the IdealGasMix::ct_gas object.

***********************************************************************/
void Reaction_set::ct_equilibrate_HP( ) const {

#ifdef _CANTERA_VERSION

  try {
    equilibrate( *ct_gas, "HP" );
  }
  catch (CanteraError) {
    Cantera::showErrors();
  }

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_parse_mass_string


void Reaction_set::ct_equilibrate_TP( ) const {

#ifdef _CANTERA_VERSION

  try {
    equilibrate( *ct_gas, "TP" );
  }
  catch (CanteraError) {
    Cantera::showErrors();
  }

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION

} // end of ct_parse_mass_string


/***********************************************************************

  This function computes the composition of a hydrocarbon (CxHy), O2, 
  and N2 mixture.  Cantera is used to convert molar fraction to mass
  fraction.

***********************************************************************/
void Reaction_set::ct_composition( const string& fuel_species, 
				   const double &phi,
				   double* massFracs ) const {

#ifdef _CANTERA_VERSION

  // first, compute the stoichiometric fuel air ratio
  double C_atoms=ct_gas->nAtoms(ct_gas->speciesIndex(fuel_species), ct_gas->elementIndex("C"));
  double H_atoms=ct_gas->nAtoms(ct_gas->speciesIndex(fuel_species), ct_gas->elementIndex("H"));
  double ax=C_atoms+H_atoms/4.0;
  double fa_stoic=1.0/(4.76*ax);

  // determine the composition
  int nsp = ct_gas->nSpecies();  
  double* x = new double[nsp];
  for(int k=0; k<nsp; k++){
    if(k==ct_gas->speciesIndex(fuel_species)){ x[k]=phi; }
    else if(k==ct_gas->speciesIndex("O2")){    x[k]=1.00*ax; }
    else if(k==ct_gas->speciesIndex("N2")){    x[k]=3.76*ax; }
    else{ x[k]=0.0; }
  }

  // compute composition 
  // -> why do it yourself when you can get someone else to do it
  ct_gas->setMoleFractions(x);
  for(int k=0;k<nsp;k++) massFracs[k] = ct_gas->massFraction(k);

  // clean up memory
  delete[] x;

#else
  cout<<"\n CODE NOT COMPILED WITH CANTERA!";
  cout<<"\n YOU SHOULD NOT BE HERE!";
  exit(-1);

#endif //_CANTERA_VERSION


} // end of ct_composition

