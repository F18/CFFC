
#ifndef _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED
#include "NavierStokes3DPolytropicInput.h"
#endif // _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED

/********************************************************
 * Routine: Set_Reference_Solution_States               *
 *                                                      *
 * Assigns default values to the reference solution     *
 * states.                                              *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Set_Reference_Solution_States(void) {

    // Set gas type.
    Wo.setgas(Gas_Type);
    Uo.setgas(Gas_Type);
    
    
    
    // Set reference primitive solution states.
    Wo.rho = Pressure/(Wo.R*Temperature); 
    Wo.p = Pressure;	
    Wo.v.zero();
    
    if(i_Original_ICs == IC_NOT_DEFINED) {
        cerr << "Need to define input parameter Original_ICs_Type";
    } else {
        if(i_Original_ICs == IC_CHANNEL_FLOW) {
            double v_max;
            double *Pressure_Gradient_ptr, *Length_ptr, *Velocity_ptr;
            switch (Grid_IP.i_Grid) {
            case GRID_CHANNEL_ZDIR :
                Pressure_Gradient_ptr = &Pressure_Gradient.z;
                Length_ptr = &Grid_IP.Box_Height;
                Velocity_ptr = &Wo.v.z;
                break;
            case GRID_CHANNEL_XDIR :
                Pressure_Gradient_ptr = &Pressure_Gradient.x;
                Length_ptr = &Grid_IP.Box_Length;
                Velocity_ptr = &Wo.v.x;
                break;
            case GRID_CHANNEL_YDIR :
                Pressure_Gradient_ptr = &Pressure_Gradient.y;
                Length_ptr = &Grid_IP.Box_Width;
                Velocity_ptr = &Wo.v.y;
                break;
            }

            /* Change grid and pressure gradient to match the given Reynolds_Number and Mach_Number */                    
            if(Mach_Number != ZERO) {
                v_max = Mach_Number * Wo.a();
                Pressure_Gradient.zero();
                *Length_ptr = Reynolds_Number * Wo.nu() / v_max;
                *Velocity_ptr = v_max;
                *Pressure_Gradient_ptr = -TWO*Wo.mu()*v_max/sqr( HALF*(*Length_ptr) );
            } else { // Mach_Number == ZERO
                /* The given Pressure_Gradient is used together with the given Grid_parameters
                to calculate the Mach_Number and Reynolds_Number */
                v_max = - *Pressure_Gradient_ptr/(TWO*Wo.mu()) * sqr(0.5 * (*Length_ptr));
                *Velocity_ptr = v_max;
                Mach_Number = v_max/Wo.a();
                Reynolds_Number = v_max* (*Length_ptr)/Wo.nu();
            }
            
            
        } else if(i_Original_ICs == IC_VISCOUS_FLAT_PLATE) {
            Wo.v.x = Mach_Number * Wo.a();
            Grid_IP.Plate_Length = Reynolds_Number * Wo.nu() / Wo.v.x;
            Grid_IP.Box_Width = TWO*Grid_IP.Plate_Length;
            Grid_IP.Box_Height = HALF*Grid_IP.Plate_Length;
            Grid_IP.Box_Length = Grid_IP.Plate_Length;
        }
      
    }
                          
    
    // Set reference conservative solution state                      
    Uo = Wo.U();  
    
}

/********************************************************
 * Routine: Read_Reference_Solution_States              *
 *                                                      *
 * Read in the reference solution states.               *
 *                                                      *
 ********************************************************/
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file) {

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
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file) {

    // Write gas type.
    restart_file << Gas_Type << endl;

    // Write reference solution states.
    restart_file << Wo << endl << Uo << endl;

}
