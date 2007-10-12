/*!\file CENO_ExecutionMode.cc
  \brief Initialize the flags that control the execution of CENO high-order reconstruction. */

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "CENO_ExecutionMode.h"

short CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
short CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
short CENO_Execution_Mode::CENO_DROP_ORDER = ON; // produce monotone solutions
short CENO_Execution_Mode::CENO_Padding = OFF; // no padding
short CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
short CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights


//! Set all flags to default values
// add all flag default values to this function
void CENO_Execution_Mode::SetDefaults(void){
  
  USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
  CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
  CENO_DROP_ORDER = ON; // produce monotone solutions
  CENO_Padding = OFF; // no padding
  CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
  CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights
}
