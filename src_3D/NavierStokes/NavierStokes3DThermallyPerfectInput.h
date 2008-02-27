/* Navierstokes3DThermallyPerfectInput.h:  Header file defining various 
                                           specializations for 
                                           Navierstokes3d_ThermallyPerfect 
                                           solution input parameter class. */

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_INPUT_INCLUDED
#define _NAVIERSTOKES3D_THERMALLYPERFECT_INPUT_INCLUDED

// Include required CFFC header files

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED
#include "NavierStokes3DThermallyPerfectState.h"
#endif // NAVIERSTOKES3D_THERMALLYPERFECT_STATE_INCLUDED

/* Define the specializations. */

//! Sets values of reference solution states
template<>
void Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                      NavierStokes3D_ThermallyPerfect_cState>::Set_Reference_Solution_States(void);

//! Read in the reference solution states
template<>
void Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                      NavierStokes3D_ThermallyPerfect_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                      NavierStokes3D_ThermallyPerfect_cState>::Write_Reference_Solution_States(ostream &restart_file);

//! Input-output operators
template<>
void Input_Parameters<NavierStokes3D_ThermallyPerfect_pState, 
                      NavierStokes3D_ThermallyPerfect_cState>::Output_Solution_Type(ostream &out_file) const;

#endif // _NAVIERSTOKES3D_THERMALLYPERFECT_INPUT_INCLUDED
