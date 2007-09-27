/* TurbulenceModellingInput.cc: Definition of Turbulence_Modelling_Input_Parameters 
                                class member functions. */

/* Include the TurbulenceModellingInput header file. */

#ifndef _TURBULENCEMODEL_INPUT_INCLUDED
#include "TurbulenceModellingInput.h"
#endif // _TURBULENCEMODEL_INPUT_INCLUDED

/* Define member functions. */

/************************************************************************************
 * Turbulence_Modelling_Input_Parameters::Broadcast -- Broadcast to all processors. *
 ************************************************************************************/
void Turbulence_Modelling_Input_Parameters::Broadcast(void) {

#ifdef _MPI_VERSION
    MPI::COMM_WORLD.Bcast(&(Near_Wall_Boundary_Treatment),
                          1,
                          MPI::INT, 0);
#endif

}

/********************************************************************************************
 * Turbulence_Modelling_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input. *
 ********************************************************************************************/
int Turbulence_Modelling_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
                                                                              stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;

  if (strcmp(code, "Near_Wall_Boundary_Treatment") == 0) {
    i_command = 6001;
    value >> Near_Wall_Boundary_Treatment;
    if  (Near_Wall_Boundary_Treatment < 0 ||
        Near_Wall_Boundary_Treatment > 2) i_command = INVALID_INPUT_VALUE;

  } else {
    i_command = INVALID_INPUT_CODE;

  } /* endif */

  return i_command;
  
}

/******************************************************************************
 * Turbulence_Modelling_Input_Parameters::Check_Inputs -- Check input values. *
 ******************************************************************************/
int Turbulence_Modelling_Input_Parameters::Check_Inputs(void) {

  // Input parameters are consistent.  Exit successfully.
  return 0;

}

/***************************************************************************
 * Turbulence_Modelling_Input_Parameters -- Input-output operators.        *
 ***************************************************************************/
ostream &operator << (ostream &out_file,
                      const Turbulence_Modelling_Input_Parameters &IP) {

  IP.Output(out_file);
  return (out_file);

}

istream &operator >> (istream &in_file,
                      Turbulence_Modelling_Input_Parameters &IP) {

  return in_file;

}

void Turbulence_Modelling_Input_Parameters::Output(ostream &out_file) const {

  if (Near_Wall_Boundary_Treatment == 0) {
     out_file << "\n  -> Near Wall Turbulent BC Treatment: Automatic";
  } else {
     if (Near_Wall_Boundary_Treatment == 2) {
        out_file << "\n  -> Near Wall Turbulent BC Treatment: Direct Integration";
     } else {
        out_file << "\n  -> Near Wall Turbulent BC Treatment: Wall Functions";
     } /* endif */
  } /* endif */

}
