/****************** Reactions.cc ************************************
  This class defines the Reaction class and for the 
   Multiple Species solution.

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
  
  /****************** METHANE +  HYDROGEN *******************************/   
  // 2 step CH4 mechanism +  2step H2 and air mechanism (Masri bluff body fuel CH4/H2) 
  else if(react == "CH4H2_2STEP"){
     //flag for hardcoded
     reactset_flag =CH4H2_2STEP;
     //set number of reactions & species
     num_reactions = 4;
     num_species = 8;
     num_react_species = 7;
     //allocate memory 
     reactions = new React_data[num_reactions];
     //set reaction coefficients
     reactions[0].set_data("CH4_2step_1",2.6*1.443e13,0.0,(48400.0*CAL_TO_JOULE),0);
     reactions[1].set_data("CH4_2step_2",pow(10.0,14.6),5.0e8,0.0,(40000.0*CAL_TO_JOULE),0);
     reactions[2].set_data("H2O2_2step_1",1.0,-10.0,4865.0*CAL_TO_JOULE,0);
     reactions[3].set_data("H2O2_2step_2",1.0,-13.0,42000.0*CAL_TO_JOULE,1);
     
     //setup species list 
     species = new string[num_species];
     species[0] = "CH4";
     species[1] = "O2";
     species[2] = "CO2";
     species[3] = "H2O";
     species[4] = "CO"; 
     species[5] = "H2";
     species[6] = "OH";
     species[7] = "N2";

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



