/* Euler3DPolytropicInput.h:  Header file defining various specializations for 
                              Euler3D_Polytropic solution input parameter class. */

#ifndef _EULER3D_POLYTROPIC_INPUT_INCLUDED
#define _EULER3D_POLYTROPIC_INPUT_INCLUDED

/* Include related CFFC header files. */

#ifndef _INPUT_INCLUDED
#include "../CFD/Input.h"
#endif // INPUT_INCLUDED

#ifndef _EULER3D_POLYTROPIC_STATE_INCLUDED
#include "Euler3DPolytropicState.h"
#endif // EULER3D_POLYTROPIC_STATE_INCLUDED

/* Define the specializations. */


//! Deallocate static data of reference solution states
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Deallocate_Static(void);


//! Sets values of reference solution states
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Set_Reference_Solution_States(void);

//! Read in the reference solution states
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<Euler3D_Polytropic_pState, 
                      Euler3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file);

#endif // _EULER3D_POLYTROPIC_INPUT_INCLUDED
