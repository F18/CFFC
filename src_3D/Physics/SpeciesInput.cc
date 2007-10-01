/* SpeciesInput.cc: Definition of Species_Input_Parameters class member functions. */

/* Include the SpeciesInput header file. */

#ifndef _SPECIESINPUT_INCLUDED
#include "SpeciesInput.h"
#endif // _SPECIESINPUT_INCLUDED

/* Define member functions. */

/*************************************************************
 * Species_Input_Parameters:Allocate -- Create Memory.       *
 *************************************************************/
void Species_Input_Parameters::Allocate(const int Ns) {

   if (Ns >= 1 && !Allocated) {
      num_species = Ns;

      Multispecies = new char*[num_species];
      for (int i = 0; i < num_species; i++) {
         Multispecies[i] = new char[SPECIES_INPUT_PARAMETER_LENGTH];
      } /* endfor */
      multispecies = new string[num_species]; 
      mass_fractions = new double[num_species];
      Schmidt = new double[num_species];

      for (int i = 0; i < num_species; i++) {
        mass_fractions[i] = ZERO;
        Schmidt[i] = Global_Schmidt;
      } /* endfor */
      mass_fractions[0] = ONE;

      Allocated = 1;
   } /* endif */

}

/*************************************************************
 * Species_Input_Parameters:Allocate -- Delete Memory.       *
 *************************************************************/
void Species_Input_Parameters::Deallocate(void) {

   if (num_species >= 1 && Allocated) {
      for (int i = 0; i < num_species; i++) {
         delete[] Multispecies[i]; Multispecies[i]=NULL;
      } /* endfor */
      delete[] Multispecies; Multispecies=NULL;
      delete[] multispecies; multispecies=NULL;
      delete[] mass_fractions; mass_fractions=NULL;
      delete[] Schmidt; Schmidt = NULL;
   
      num_species = 0;

      Allocated = 0;
   } /* endif */

}

/***************************************************************************
 * Species_Input_Parameters::Broadcast -- Broadcast to all processors.     *
 ***************************************************************************/
void Species_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
   // Reaction Name:
   MPI::COMM_WORLD.Bcast(Reaction_Mechanism_Name,
                         SPECIES_INPUT_PARAMETER_LENGTH,
                         MPI::CHAR, 0);

   // Cantera Reaction Mechanism:
   MPI::COMM_WORLD.Bcast(Cantera_Mech_Name,
                         SPECIES_INPUT_PARAMETER_LENGTH,
                         MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(Cantera_Mech_File,
                         SPECIES_INPUT_PARAMETER_LENGTH,
                         MPI::CHAR, 0);

   // Delete current dynamic memory before changing num_species.
   if (!CFFC_Primary_MPI_Processor()) {
      Deallocate();
   } /* endif */

   // Number of Species:
   int ns;
   if (CFFC_Primary_MPI_Processor()) {
      ns=num_species;
   } /* endif */
   MPI::COMM_WORLD.Bcast(&(ns),
                         1,
                         MPI::INT, 0);

   // Allocate required dynamic memory for species data.
   if(!CFFC_Primary_MPI_Processor()) {
      Allocate(ns);
   } /* endif */

   // Species names & mass fractions:
   for (int i =0; i < num_species; i++) {
      MPI::COMM_WORLD.Bcast(Multispecies[i],
			    SPECIES_INPUT_PARAMETER_LENGTH,
			    MPI::CHAR, 0);
      MPI::COMM_WORLD.Bcast(&(mass_fractions[i]),
			    1,
			    MPI::DOUBLE, 0);
      MPI::COMM_WORLD.Bcast(&(Schmidt[i]),
			    1,
			    MPI::DOUBLE, 0);
   } /* endfor*/

   // Set additional reaction and species parameters.
   if (!CFFC_Primary_MPI_Processor()) {
      reaction_mechanism_name = Reaction_Mechanism_Name;
      cantera_mech_name = Cantera_Mech_Name;
      cantera_mech_file = Cantera_Mech_File;
      for (int i = 0; i < num_species; i++) {
         multispecies[i] = Multispecies[i];
      } /* endfor */
   } /* endif */

   // Global Schmidt number:
   MPI::COMM_WORLD.Bcast(&(Global_Schmidt),
                         1,
                         MPI::INT, 0);

   // Transport data:
   MPI::COMM_WORLD.Bcast(Transport_Data_Type, 
        		 SPECIES_INPUT_PARAMETER_LENGTH, 
			 MPI::CHAR, 0);
   MPI::COMM_WORLD.Bcast(&(i_transport_data_type), 
			 1, 
			 MPI::INT, 0);
   if (!CFFC_Primary_MPI_Processor()) {
     transport_data_type = Transport_Data_Type;
   } /* endif */
#endif

}

/*******************************************************************************
 * Species_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 *******************************************************************************/
int Species_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                                 stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  Reaction_set reaction_mechanism;

  /*************************************/
  /**** REACTIONS SET FOR HARDCODED ****/
  /*************************************/
  if (strcmp(code, "Reaction_Mechanism") == 0) {
    i_command = 5001;
    value >> reaction_mechanism_name;
    strcpy(Reaction_Mechanism_Name, reaction_mechanism_name.c_str());

    Deallocate();
    reaction_mechanism.set_reactions(reaction_mechanism_name);
    num_species = reaction_mechanism.num_species;       
    Allocate(reaction_mechanism.num_species);
       
    // Get species and load appropriate data
    for (int i = 0; i < num_species; i++) {
       multispecies[i] = reaction_mechanism.species[i];
       strcpy(Multispecies[i], multispecies[i].c_str());
       Schmidt[i] = Global_Schmidt;
    } /* endfor */
 
  /*************************************/
  /***** REACTIONS SET FOR CANTERA *****/
  /*************************************/
  } else if (strcmp(code, "Cantera_Reaction_Mechanism_Name") == 0) {
    i_command = 5002;
    reaction_mechanism_name = "CANTERA";
    strcpy(Reaction_Mechanism_Name, reaction_mechanism_name.c_str());
    value >> cantera_mech_name;
    strcpy(Cantera_Mech_Name, cantera_mech_name.c_str());

    // Using existing mechanism file , load the reaction mechanism.
    if (cantera_mech_file != "none") {
       Deallocate();
       reaction_mechanism.ct_load_mechanism(cantera_mech_file, cantera_mech_name);   
       num_species = reaction_mechanism.num_species;
       Allocate(reaction_mechanism.num_species);

       // Get species and load appropriate data
       for (int i = 0; i < num_species; i++) {
          multispecies[i] = reaction_mechanism.species[i];
          strcpy(Multispecies[i], multispecies[i].c_str());
          Schmidt[i] = Global_Schmidt;
       } /* endfor */
    } /* endif */

  } else if (strcmp(code, "Cantera_Reaction_Mechanism_File") == 0) {
    i_command = 5003;
    reaction_mechanism_name = "CANTERA";
    strcpy(Reaction_Mechanism_Name, reaction_mechanism_name.c_str());
    value >> cantera_mech_file;
    strcpy(Cantera_Mech_File, cantera_mech_file.c_str());

    // Using existing mechanism name, load the reaction mechanism.
    if (cantera_mech_name != "none") {
       Deallocate();
       reaction_mechanism.ct_load_mechanism(cantera_mech_file, cantera_mech_name);   
       num_species = reaction_mechanism.num_species;
       Allocate(reaction_mechanism.num_species);

       // Get species and load appropriate data
       for (int i = 0; i < num_species; i++) {
          multispecies[i] = reaction_mechanism.species[i];
          strcpy(Multispecies[i], multispecies[i].c_str());
          Schmidt[i] = Global_Schmidt;
       } /* endfor */
    } /* endif */

  /******************************************/
  /**** NON REACTING, BUT MULTIPLE GASES ****/
  /******************************************/
  } else if (strcmp(code, "Species") == 0) { 
    i_command = 5004;
    reaction_mechanism_name ="NO_REACTIONS";
    strcpy(Reaction_Mechanism_Name, reaction_mechanism_name.c_str());

    Deallocate();
    value >> num_species;   
    Allocate(num_species);
       
    // Obtain species names.
    for (int i = 0; i < num_species; i++) {
       value >> multispecies[i];
    } /* endfor */

    // Set appropriate species data
    for (int i = 0; i < num_species; i++){
       strcpy(Multispecies[i], multispecies[i].c_str());
       Schmidt[i] = Global_Schmidt;
    } /* endfor */
 
  /***************************************/
  /**** REACTIONS FOR USER DEFINED *******/ 
  /***************************************/
  } else if (strcmp(code, "User_Reaction_Mechanism") == 0) { 
    i_command = 5005;
    cout << endl << code 
         << "\n not currently available in freeware version :)\n";
    i_command = INVALID_INPUT_VALUE;

  /****************************/
  /**** Schmidt Numbers *******/ 
  /****************************/
  } else if (strcmp(code, "Schmidt") == 0) {
    i_command = 5006; 
    value >> Global_Schmidt;
    if (Global_Schmidt < 0) {
       i_command = INVALID_INPUT_VALUE;
    } else {
       for (int i = 0; i < num_species; i++){
          Schmidt[i] = Global_Schmidt;
       } /* endfor */
    } /* endif */

  } else if (strcmp(code, "Schmidt_Numbers") == 0) {
    i_command = 5007; 
    for (int i=0; i < num_species; i++){
       value >> Schmidt[i];	 	 
    } /* endfor */
 
  /***************************/
  /**** Mass Fractions *******/ 
  /***************************/
  } else if (strcmp(code, "Mass_Fractions") == 0) {
    i_command = 5008; 
    if (reaction_mechanism_name != "CANTERA") {
       // Get Initial Mass Fractions from user 
       double temp = ZERO;
       for (int i = 0; i < num_species; i++){
          value >> mass_fractions[i];
          temp += mass_fractions[i];
       } /*endfor */
       //check to make sure mass fractions adds to 1
       if (temp < ONE-MICRO || temp > ONE+MICRO) { 
          cout << "\n Error: input mass fractions have a sum of " 
               << temp
               << ". Should sum to 1.";
          cout << "\n mass_fractions: ";
          for (int i = 0; i < num_species; i++){
	     cout << mass_fractions[i] << " ";
          } /*endfor */
          i_command = INVALID_INPUT_VALUE;
       } /* endif */
    } else {
       // Here we use Cantera to parse the string of the form:
       //       CH4:0.5, O2:0.5
       // All other species will be assumed to have 0 mass fractions.  
       // Cantera also normalizes the mass fractions to sum to unity.  
       // Returns them in an array.
       string mass_fractions_input;
       char Mass_Fractions_Input[SPECIES_INPUT_PARAMETER_LENGTH];
       mass_fractions_input = value.str();
       strcpy(Mass_Fractions_Input, mass_fractions_input.c_str());
       reaction_mechanism.ct_load_mechanism(cantera_mech_file, cantera_mech_name);
       reaction_mechanism.ct_parse_mass_string(Mass_Fractions_Input, mass_fractions);
    } /* endif */

  /********* TRANSPORT DATA **********************************************
  // Transport-Lennard-Jones  - Use Lennard-Jones potentials or 
  // Transport-NASA - use Nasa curvefits
  // If not specified, default = NASA
  ************************************************************************/
  } else if (strcmp(code, "Transport_Data_Type") == 0) {
    i_command = 5009; 
    value >> transport_data_type;
    strcpy(Transport_Data_Type, transport_data_type.c_str());
    if (transport_data_type == "Transport-NASA") {
       i_transport_data_type = TRANSPORT_NASA;
    } else if (transport_data_type == "Transport-Lennard-Jones") {
       i_transport_data_type = TRANSPORT_LENNARD_JONES;
    } else {
       i_command = INVALID_INPUT_VALUE;
    } // end if

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/***************************************************************************
 * Species_Input_Parameters::Check_Inputs -- Check input values.           *
 ***************************************************************************/
int Species_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * Species_Input_Parameters -- Input-output operators.                     *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Species_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      Species_Input_Parameters &IP) {

  return in_file;

}

void Species_Input_Parameters::Output(ostream &out_file) const { 

  out_file << "\n  -> Reaction Mechanism: " 
           << reaction_mechanism_name;
  if (reaction_mechanism_name == "CANTERA") {
     out_file << "\n  -> Cantera Reaction Mechanism Name: " 
              << cantera_mech_name;
     out_file << "\n  -> Cantera Reaction Mechanism File: " 
              << cantera_mech_file;
  } /* endif */
  out_file << "\n  -> Number of Species: " 
           << num_species;
  out_file << "\n  -> Mixture components: ";
  for (int i = 0; i < num_species; i++) {
     out_file << multispecies[i];
     if (i != num_species-1) out_file << ", ";
  } /* endfor */
  out_file << "\n  -> Initial mass fractions: ";
  for (int i = 0; i < num_species; i++) {
     out_file  << "c[" << multispecies[i] << "]= ";
     out_file  << mass_fractions[i];
     if (i != num_species-1) out_file << ", ";
  } /* endfor */
  out_file << "\n  -> Species Schmidt Numbers: ";
  for (int i = 0; i < num_species; i++) {
     out_file << "Sc[" << multispecies[i] << "]= ";
     out_file << Schmidt[i];
     if (i != num_species-1) out_file << ", ";
  } /* endfor */
  out_file << "\n  -> Transport Data Type: "
           << transport_data_type;

}
