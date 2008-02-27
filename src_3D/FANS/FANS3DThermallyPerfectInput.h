/* FANS3DThermallyPerfectInput.h:  Header file defining various specializations for 
                                   FANS3D_ThermallyPerfect solution input parameter class. */

#ifndef _FANS3D_THERMALLYPERFECT_INPUT_INCLUDED
#define _FANS3D_THERMALLYPERFECT_INPUT_INCLUDED

/* Include related CFFC header files. */

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _FANS3D_THERMALLYPERFECT_STATE_INCLUDED
#include "FANS3DThermallyPerfectState.h"
#endif // _FANS3D_THERMALLYPERFECT_STATE_INCLUDED

/* Define the specializations. */

//! Sets values of reference solution states
template<>
void Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                      FANS3D_ThermallyPerfect_KOmega_cState>::Set_Reference_Solution_States(void);

//! Read in the reference solution states
template<>
void Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                      FANS3D_ThermallyPerfect_KOmega_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                      FANS3D_ThermallyPerfect_KOmega_cState>::Write_Reference_Solution_States(ostream &restart_file);

//! Input-output operators
template<>
void Input_Parameters<FANS3D_ThermallyPerfect_KOmega_pState, 
                      FANS3D_ThermallyPerfect_KOmega_cState>::Output_Solution_Type(ostream &out_file) const;

#endif // _FANS3D_THERMALLYPERFECT_INPUT_INCLUDED
