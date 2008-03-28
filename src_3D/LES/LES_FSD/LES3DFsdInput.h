/* LES3DFsdInput.h:  Header file defining various specializations for 
                     LES3DFsd solution input parameter class. */

#ifndef _LES3DFSD_INPUT_INCLUDED
#define _LES3DFSD_INPUT_INCLUDED

// Include required CFFC header files

#ifndef _INPUT_INCLUDED
#include "../../CFD/Input.h"
#endif // _INPUT_INCLUDED

#ifndef _LES3DFSD_STATE_INCLUDED
#include "LES3DFsdState.h"
#endif // _LES3DFSD_STATE_INCLUDED

/* Define the specializations. */

//! Deallocates static data of reference solution states
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Deallocate_Static(void);

//! Sets values of reference solution states
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Set_Reference_Solution_States(void);
  
//! Read in the reference solution states
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Read_Reference_Solution_States(istream &restart_file);

//! Write out the reference solution states
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Write_Reference_Solution_States(ostream &restart_file);

//! Input-output operators
template<>
void Input_Parameters<LES3DFsd_pState, 
                      LES3DFsd_cState>::Output_Solution_Type(ostream &out_file) const;

#endif // _LES3DFSD_INPUT_INCLUDED


