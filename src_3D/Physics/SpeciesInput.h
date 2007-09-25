/* SpeciesInput.h:  Header file declaring multi-species 
                    gaseous mixture input class. */

#ifndef _SPECIESINPUT_INCLUDED
#define _SPECIESINPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#ifndef _SPECIES_INCLUDED
#include "Species.h"
#endif // _SPECIES_INCLUDED

#ifndef _NASARP1311_DATA_INCLUDED
#include "NASAData/NASARP1311data.h"
#endif // _NASARP1311_DATA_INCLUDED

#ifndef _REACTIONS_INCLUDED
#include "../Reactions/Reactions.h"
#endif // _REACTIONS_INCLUDED

#define SPECIES_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: Species_Input_Parameters
 *
 * @brief Input Parameters for multi-species gaseous mixtures.
 *
 * This class defines and handles the input variables related to the
 * treatment of mulit-species gaseous mixtures.
 *
 */
class Species_Input_Parameters{
 private:
 public:
  //@{ @name Reaction data input parameters:
  //! Name of reaction mechanism.
  string reaction_mechanism_name;
  char Reaction_Mechanism_Name[SPECIES_INPUT_PARAMETER_LENGTH];
  //@}

  //@{ @name Cantera input parameters:
  //! Mechanism name:
  string cantera_mech_name;
  char Cantera_Mech_Name[SPECIES_INPUT_PARAMETER_LENGTH];
  //! Mechanism file:
  string cantera_mech_file;
  char Cantera_Mech_File[SPECIES_INPUT_PARAMETER_LENGTH];
  //@}

  //@{ @name Speceis input parameters:
  //! Number of species in mixture:
  int num_species;
  //! Names of each species in mixture:
  char **Multispecies;
  string *multispecies;
  //! Mass fractions of each species in mixture:
  double *mass_fractions;
  //! Schmidt numbers of each species in mixture:
  double *Schmidt;
  double Global_Schmidt;  //depricated, use each individual Schmidt's
  //@}

  //@{ @name Transport data type:
  string transport_data_type;
  char Transport_Data_Type[SPECIES_INPUT_PARAMETER_LENGTH];
  int  i_transport_data_type; 
  //@}

  //@{ @name Variable allocation indicator for species input parameters:
  int Allocated;
  //@}

  //@{ @name Constructor and desctructor
  //! Constructor (assign default values):
  Species_Input_Parameters() {
     Allocated = 0;
     //Default: multi-species with no reactions.
     reaction_mechanism_name ="NO_REACTIONS";
     strcpy(Reaction_Mechanism_Name, reaction_mechanism_name.c_str());
     cantera_mech_name = "none";
     cantera_mech_file = "none";
     strcpy(Cantera_Mech_Name, cantera_mech_name.c_str());
     strcpy(Cantera_Mech_File, cantera_mech_file.c_str());
     //Use air with 79% N2, and 21% 02 by volume.
     Allocate(2);
     multispecies[0] = "N2"; 
     multispecies[1] = "O2"; 
     strcpy(Multispecies[0], multispecies[0].c_str());
     strcpy(Multispecies[1], multispecies[1].c_str());
     mass_fractions[0] = 0.765; 
     mass_fractions[1] = 0.235;
     Global_Schmidt = ONE;
     Schmidt[0] = Global_Schmidt;
     Schmidt[1] = Global_Schmidt;
     //Set transport data type.
     transport_data_type = "Transport-NASA";
     strcpy(Transport_Data_Type, transport_data_type.c_str());
     i_transport_data_type = TRANSPORT_NASA;
  }

  Species_Input_Parameters(const int Ns){
     Allocate(Ns);    
  }

  //! Destructor:
  ~Species_Input_Parameters(void) {
     Deallocate();
  }
  //@}

  //@{ @name Other Member functions.
  //! Allocate required memory:
  void Allocate(const int Ns);
  //! Deallocate require memory:
  void Deallocate(void);
  //! Broadcast input parameters to all processors:
  void Broadcast(void);
  //! Parse next input line:
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
  //! Check validity of specified input parameters:
  int Check_Inputs(void);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
		               const Species_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Species_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  //@}

};

#endif // _SPECIESINPUT_INCLUDED
