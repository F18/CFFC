/* Euler3DThermallyPerfectInput.h:  Header file defining various specializations for 
                                    Euler3D_ThermallyPerfect solution input parameter class. */

#ifndef _EULER3D_THERMALLYPERFECT_INPUT_INCLUDED
#define _EULER3D_THERMALLYPERFECT_INPUT_INCLUDED

/* Include related CFFC header files. */

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _EULER3D_THERMALLYPERFECT_STATE_INCLUDED
#include "Euler3DThermallyPerfectState.h"
#endif // _EULER3D_THERMALLYPERFECT_STATE_INCLUDED

/* Define the specializations. */

/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Set_Reference_Solution_States(void) {

    // Load reaction mechanism data.
    if (Species_IP.reaction_mechanism_name != "CANTERA") {
       Wo.React.set_reactions(Species_IP.reaction_mechanism_name);
    } else {
       Wo.React.ct_load_mechanism(Species_IP.cantera_mech_file, 
                                  Species_IP.cantera_mech_name);
    } /* endif */
      
    // Set species types if non-reacting mixture.
    if (Wo.React.reactset_flag == NO_REACTIONS) {
       Wo.React.set_species(Species_IP.multispecies, 
                            Species_IP.num_species);
    } /* endif */

    // Obtain and set physical/transport data for component gases
    // as well as set default default initial values.
    Wo.set_species_data(Species_IP.num_species,
                        Species_IP.multispecies,
                        CFFC_Path,
	                Debug_Level,
                        Mach_Number_Reference,
                        Species_IP.Schmidt,
                        Species_IP.i_transport_data_type); 
    Uo.set_species_data(Species_IP.num_species,
                        Species_IP.multispecies,
                        CFFC_Path,
                        Debug_Level,
                        Mach_Number_Reference,
                        Species_IP.Schmidt,
                        Species_IP.i_transport_data_type);
    Wo.set_initial_values(Species_IP.mass_fractions);
    Uo.set_initial_values(Species_IP.mass_fractions);

    // Set reference solution states.
    Wo.rho = Pressure/(Wo.Rtot()*Temperature); 
    Wo.p = Pressure;	
    Wo.v.zero();
    Wo.v.x = Mach_Number*Wo.a()*cos(TWO*PI*Flow_Angle/360.00);
    Wo.v.y = Mach_Number*Wo.a()*sin(TWO*PI*Flow_Angle/360.00);
    Uo = Wo.U();  

}

/********************************************************
 * Routine: Read_Reference_Solution_States              *
 *                                                      *
 * Read in the reference solution states.               *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Read_Reference_Solution_States(istream &restart_file) {

    char line[256];

    // Read name of reaction mechanism.
    restart_file.getline(line,sizeof(line)); 
    restart_file >> Species_IP.reaction_mechanism_name;
    restart_file.getline(line,sizeof(line)); 
    restart_file >> Species_IP.cantera_mech_name;
    restart_file.getline(line,sizeof(line)); 
    restart_file >> Species_IP.cantera_mech_file;

    strcpy(Species_IP.Reaction_Mechanism_Name, Species_IP.reaction_mechanism_name.c_str());
    strcpy(Species_IP.Cantera_Mech_Name, Species_IP.cantera_mech_name.c_str());
    strcpy(Species_IP.Cantera_Mech_File, Species_IP.cantera_mech_file.c_str());

    // Set skipping of white spaces.       
    restart_file.setf(ios::skipws);

    // Read the number of species.
    Species_IP.Deallocate();
    restart_file >> Species_IP.num_species; 
    Species_IP.Allocate(Species_IP.num_species);
    
    // Read in species data.
    for (int i = 0; i < Species_IP.num_species; i++) {
       restart_file >> Species_IP.multispecies[i];
       strcpy(Species_IP.Multispecies[i], Species_IP.multispecies[i].c_str());
    } /* endfor */ 

    // Read species initial mass fractions and Schmidt numbers.
    for (int i = 0; i < Species_IP.num_species; i++) {
       restart_file >> Species_IP.mass_fractions[i] >> Species_IP.Schmidt[i];
    } /* endfor */

    // Stop skipping white spaces.
    restart_file.unsetf(ios::skipws);

    // Set reference solution states.
    Set_Reference_Solution_States();

    // Read in reference solution states.
    restart_file >> Wo >> Uo;

}

/********************************************************
 * Routine: Write_Reference_Solution_States             *
 *                                                      *
 * Write out the reference solution states.             *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Write_Reference_Solution_States(ostream &restart_file) {

    // Write name of reaction mechanism.
    restart_file << Species_IP.reaction_mechanism_name << endl;
    restart_file << Species_IP.cantera_mech_name << endl;
    restart_file << Species_IP.cantera_mech_file << endl;

    // Write out species data.
    restart_file << Species_IP.num_species << " ";
    for (int i = 0; i < Species_IP.num_species; i++){ 
       restart_file << Species_IP.multispecies[i] << " ";
    } /* endfor */
    restart_file << endl;

    // Write species initial mass fractions and Schmidt numbers.
    for (int i = 0; i < Species_IP.num_species; i++) {
      restart_file << Species_IP.mass_fractions[i] << " " 
                   << Species_IP.Schmidt[i] << endl;
    } /* endfor */

    // Write reference solution states.
    restart_file << Wo << endl << Uo << endl;

}

/******************************************************************************
 * Euler3D_ThermallyPerfect_Input_Parameters -- Input-output operators.       *
 ******************************************************************************/
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Output_Solution_Type(ostream &out_file) const {

    // Multi-species Input Parameters:
    out_file << "\n\n Thermally Perfect Reactive/Non-Reactive Gaseous Mixture";
    out_file << Species_IP;

}

#endif // _EULER3D_THERMALLYPERFECT_INPUT_INCLUDED
