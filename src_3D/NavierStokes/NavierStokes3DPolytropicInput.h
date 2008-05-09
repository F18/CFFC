/* Navierstokes3DPolytropicInput.h:  Header file defining various 
                                           specializations for 
                                           Navierstokes3D_Polytropic 
                                           solution input parameter class. */

#ifndef _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED
#define _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED

// Include required CFFC header files

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED
#include "NavierStokes3DPolytropicState.h"
#endif // NAVIERSTOKES3D_POLYTROPIC_STATE_INCLUDED

/* Define the specializations. */

//! Deallocate static data of reference solution states
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Deallocate_Static(void);

//!Sets values of reference solution states
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Set_Reference_Solution_States(void);

//! Read in the reference solution states
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<NavierStokes3D_Polytropic_pState, 
                      NavierStokes3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file);

#endif // _NAVIERSTOKES3D_POLYTROPIC_INPUT_INCLUDED
