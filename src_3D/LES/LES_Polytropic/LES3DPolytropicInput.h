/* LES3DPolytropicInput.h:  Header file defining various specializations for 
                     LES3D solution input parameter class. */

#ifndef _LES3D_POLYTROPIC_INPUT_INCLUDED
#define _LES3D_POLYTROPIC_INPUT_INCLUDED

// Include required CFFC header files

#ifndef _INPUT_INCLUDED
#include "../../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _LES3D_POLYTROPIC_STATE_INCLUDED
#include "LES3DPolytropicState.h"
#endif // _LES3D_POLYTROPIC_STATE_INCLUDED

/* Define the specializations. */


//! Deallocate static data of reference solution states
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Deallocate_Static(void);

//!Sets values of reference solution states
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Set_Reference_Solution_States(void);
  
//! Read in the reference solution states
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Write_Reference_Solution_States(ostream &restart_file);

//! Input-output operators
template<>
void Input_Parameters<LES3D_Polytropic_pState, 
                      LES3D_Polytropic_cState>::Output_Solution_Type(ostream &out_file) const;

#endif // _LES3D_INPUT_INCLUDED


