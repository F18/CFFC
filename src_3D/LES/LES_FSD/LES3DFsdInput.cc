
#ifndef _LES3DFSD_INPUT_INCLUDED
#include "LES3DFsdInput.h"
#endif // _LES3DFSD_INPUT_INCLUDED


/********************************************************
 * Routine: Deallocate_Static                           *
 *                                                      *
 * Deallocate static data of the reference solution     * 
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Deallocate_Static(void) {

  Wo.Deallocate_static(); 
  Uo.Deallocate_static();

}

/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Set_Reference_Solution_States(void) {
  
    // Load reaction mechanism data.
    if (Species_IP.reaction_mechanism_name != "CANTERA") {
       Wo.React.set_reactions(Species_IP.reaction_mechanism_name);
    } else {
       Wo.React.ct_load_mechanism(Species_IP.cantera_mech_file,
                                  Species_IP.cantera_mech_name);
    } /* endif */
      
    // Set species types if non-reacting.
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

    // Set some modelling parameters.
    Wo.set_modelling_data(Turbulence_IP.Fuel_Equivalence_Ratio,
			  Turbulence_IP.Unburnt_Fuel_Mass_Fraction,
			  Turbulence_IP.Reactants_Density,
                          Turbulence_IP.Laminar_Flame_Speed,
			  Turbulence_IP.Laminar_Flame_Thickness,
			  Turbulence_IP.Adiabatic_Flame_Temperature,
			  Turbulence_IP.Filter_Width);
    Uo.set_modelling_data(Turbulence_IP.Fuel_Equivalence_Ratio,
			  Turbulence_IP.Unburnt_Fuel_Mass_Fraction,
			  Turbulence_IP.Reactants_Density,
                          Turbulence_IP.Laminar_Flame_Speed,
			  Turbulence_IP.Laminar_Flame_Thickness,
			  Turbulence_IP.Adiabatic_Flame_Temperature,
			  Turbulence_IP.Filter_Width);
//     Uo.set_modelling_data(Turbulence_IP.Laminar_Flame_Speed,
// 			  Turbulence_IP.Laminar_Flame_Thickness,
// 			  Turbulence_IP.Thickening_Factor,
// 			  Turbulence_IP.Filter_Width);  

    // Set reference solution states.
    Wo.rho = 1.13; 
    Wo.p = Pressure;	
    Wo.v.zero();
    Wo.v.x = 0.3837;
    Wo.k = 0.0;
    Wo.C = 0.0;
    Wo.Fsd = 0.0;
    Uo = Wo.U();

}

/********************************************************
 * Routine: Read_Reference_Solution_States              *
 *                                                      *
 * Read in the reference solution states.               *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Read_Reference_Solution_States(istream &restart_file) {

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
       
    // Read in species data.
    restart_file.setf(ios::skipws);
    // Read the number of species.
    Species_IP.Deallocate();
    restart_file >> Species_IP.num_species; 
    Species_IP.Allocate(Species_IP.num_species);
    
    for (int i = 0; i < Species_IP.num_species; i++) {
       restart_file >> Species_IP.multispecies[i];
       strcpy(Species_IP.Multispecies[i], Species_IP.multispecies[i].c_str());
    } /* endfor */ 

    // Read species initial mass fractions and Schmidt numbers.
    for (int i = 0; i < Species_IP.num_species; i++) {
       restart_file >> Species_IP.mass_fractions[i] >> Species_IP.Schmidt[i];
    } /* endfor */

    // Read in modelling data.
    restart_file >> Turbulence_IP.Fuel_Equivalence_Ratio
		 >> Turbulence_IP.Unburnt_Fuel_Mass_Fraction
		 >> Turbulence_IP.Reactants_Density
		 >> Turbulence_IP.Laminar_Flame_Speed
		 >> Turbulence_IP.Laminar_Flame_Thickness
		 >> Turbulence_IP.Adiabatic_Flame_Temperature
		 >> Turbulence_IP.Filter_Width;

//     restart_file >> Turbulence_IP.Laminar_Flame_Speed 
// 		 >> Turbulence_IP.Laminar_Flame_Thickness
// 		 >> Turbulence_IP.Thickening_Factor
// 		 >> Turbulence_IP.Filter_Width;

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
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Write_Reference_Solution_States(ostream &restart_file) {

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

    // Write out modelling data.
    restart_file << Turbulence_IP.Fuel_Equivalence_Ratio      << " "
		 << Turbulence_IP.Unburnt_Fuel_Mass_Fraction  << " "
		 << Turbulence_IP.Reactants_Density           << " "
		 << Turbulence_IP.Laminar_Flame_Speed         << " "
		 << Turbulence_IP.Laminar_Flame_Thickness     << " "
		 << Turbulence_IP.Adiabatic_Flame_Temperature << " "
		 << Turbulence_IP.Filter_Width                << endl;

//     restart_file << Turbulence_IP.Laminar_Flame_Speed     << " "
// 		 << Turbulence_IP.Laminar_Flame_Thickness << " "
// 		 << Turbulence_IP.Thickening_Factor       << " "
// 		 << Turbulence_IP.Filter_Width            << endl;

    // Write reference solution states.
    restart_file << Wo << endl << Uo << endl;

}

/*************************************************************************************
 * LES3DFsd_Input_Parameters -- Input-output operators.                              *
 *************************************************************************************/
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Output_Solution_Type(ostream &out_file) const {

    // Multi-species Input Parameters:
    out_file << "\n\n Thermally Perfect Reactive/Non-Reactive Gaseous Mixture";
    out_file << Species_IP;

}
