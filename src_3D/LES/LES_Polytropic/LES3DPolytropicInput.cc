#ifndef _LES3D_POLYTROPIC_INPUT_INCLUDED
#include "LES3DPolytropicInput.h"
#endif // _LES3D_POLYTROPIC_INPUT_INCLUDED
    

/********************************************************
 * Routine: Deallocate_Static                           *
 *                                                      *
 * Deallocate static data of the reference solution     * 
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Deallocate_Static(void) {
        
}


/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3D_Polytropic_pState,LES3D_Polytropic_cState>::
Set_Reference_Solution_States(void) {
  
    // Set gas
    Wo.setgas(Gas_Type);
    Uo.setgas(Gas_Type);
    
    // Set reference solution states.
    Wo.rho = Pressure/(Wo.R*Temperature); 
    Wo.p = Pressure;	
    Wo.v.zero();
    Uo = Wo.U();
    
    // create filter
    Filter_Parameters filter;
    if (Turbulence_IP.SFS_FGR == DEFAULT){
        filter.FGR = 0.75 * ExplicitFilters_IP.FGR[Explicit_Filter_Constants::PRIMARY_FILTER];
        // Smagorinsky FGR width must be smaller than explicit FGR since it cuts off at lower wave numbers
    }
    else {
        filter.FGR = Turbulence_IP.SFS_FGR;
    }
    filter.type = ExplicitFilters_IP.Filter_Type[Explicit_Filter_Constants::PRIMARY_FILTER];
                          
    SFS_model_Parameters SFS_model;
    SFS_model.model = Turbulence_IP.i_SFS_model;
    SFS_model.smagorinsky_coefficient = Turbulence_IP.smagorinsky_coefficient;    
                          
    // Set LES parameters
    Wo.Set_LES_parameters(SFS_model,filter);
    Uo.Set_LES_parameters(SFS_model,filter);
    
}

/********************************************************
 * Routine: Read_Reference_Solution_States              *
 *                                                      *
 * Read in the reference solution states.               *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file) {

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
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file) {
    // Write gas type.
    restart_file << Gas_Type << endl;

    // Write reference solution states.
    restart_file << Wo << endl << Uo << endl;

}

/*************************************************************************************
 * LES3D_Input_Parameters -- Input-output operators.                              *
 *************************************************************************************/
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Output_Solution_Type(ostream &out_file) const {
    // Turbulence Modelling Input Parameters:
    out_file << "\n\n Turbulence Modelling Input Parameters";
    out_file << Turbulence_IP;

}
