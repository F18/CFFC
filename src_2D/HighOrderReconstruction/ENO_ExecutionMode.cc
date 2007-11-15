/*!\file ENO_ExecutionMode.cc
  \brief Initialize the flags that control the execution of ENO high-order reconstruction. */

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "ENO_ExecutionMode.h"

short ENO_Execution_Mode::USE_ENO_ALGORITHM = OFF; // ENO scheme is not used
short ENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface
short ENO_Execution_Mode::USE_PRIMITIVE_VARIABLES = OFF; // the primitive variables are calculated from the conserved ones


//! Set all flags to default values
// add all flag default values to this function
void ENO_Execution_Mode::SetDefaults(void){
  
  USE_ENO_ALGORITHM = OFF; // ENO scheme is not used
  FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface
  USE_PRIMITIVE_VARIABLES = OFF; // the primitive variables are calculated from the conserved ones
}
