
#ifndef _EULER3D_POLYTROPIC_INPUT_INCLUDED
#include "Euler3DPolytropicInput.h"
#endif // _EULER3D_POLYTROPIC_INPUT_INCLUDED


/********************************************************
 * Routine: Deallocate_Static                           *
 *                                                      *
 * Deallocate static data of the reference solution     * 
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Deallocate_Static(void) {
    
}



/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Set_Reference_Solution_States(void) {

    // Set gas type.
    Wo.setgas(Gas_Type);
    Uo.setgas(Gas_Type);

    // Set reference solution states.
    Wo.rho = Pressure/(Wo.R*Temperature); 
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
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file) {

    char line[256];

    // Read gas type.
    restart_file.getline(line,sizeof(line)); 
    restart_file >> Gas_Type;

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
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file) {

    // Write gas type.
    restart_file << Gas_Type << endl;

    // Write reference solution states.
    restart_file << Wo << endl << Uo << endl;

}

