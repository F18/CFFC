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

//! Deallocate static data of reference solution states
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState,
                      Euler3D_ThermallyPerfect_cState>::Deallocate_Static(void);

//! Set values of reference solution states
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Set_Reference_Solution_States(void);

//! Read in the reference solution states
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Write_Reference_Solution_States(ostream &restart_file);

//! Input-output operators
template<>
void Input_Parameters<Euler3D_ThermallyPerfect_pState, 
                      Euler3D_ThermallyPerfect_cState>::Output_Solution_Type(ostream &out_file) const;

#endif // _EULER3D_THERMALLYPERFECT_INPUT_INCLUDED
