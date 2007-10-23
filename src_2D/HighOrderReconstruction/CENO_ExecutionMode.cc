/*!\file CENO_ExecutionMode.cc
  \brief Initialize the flags that control the execution of CENO high-order reconstruction. */

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "CENO_ExecutionMode.h"

short CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
short CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
short CENO_Execution_Mode::CENO_DROP_ORDER = ON; // produce monotone solutions
short CENO_Execution_Mode::CENO_PADDING = OFF; // no padding
short CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
short CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights
short CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface


//! Set all flags to default values
// add all flag default values to this function
void CENO_Execution_Mode::SetDefaults(void){
  
  USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
  CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
  CENO_DROP_ORDER = ON; // produce monotone solutions
  CENO_PADDING = OFF; // no padding
  CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
  CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights
  FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void CENO_Execution_Mode::Print_Info(std::ostream & out_file){

  // output execution mode
  if (CENO_SPEED_EFFICIENT == ON){
    out_file << "\n     -> Execution Mode = " << "Speed Efficient";
  } else {
    out_file << "\n     -> Execution Mode = " << "Memory Efficient";
  }

  // output monotonicity mode
  if (CENO_DROP_ORDER == ON){
    out_file << "\n     -> Monotonicity Mode = " << "Yes (drop order)";
  } else {
    out_file << "\n     -> Monotonicity Mode = " << "No (don't drop order)";
  }

  // output padding mode
  if (CENO_PADDING == ON){
    out_file << "\n     -> Cell Padding = " << "Yes (flag adjacent cells too)";
  } else {
    out_file << "\n     -> Cell Padding = " << "No";
  }

  // output geom weighting type
  if (CENO_SQUARE_GEOM_WEIGHTING == ON){
    out_file << "\n     -> Geom Weighting = " << "Inverse of squared distance";
  } else {
    out_file << "\n     -> Geom Weighting = " << "Inverse distance";
  }

  // output interface behavior
  if (FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE == ON){
    out_file << "\n     -> Negative interface solutions = " << "Force with PWC";
  } else {
    out_file << "\n     -> Negative interface solutions = " << "Exit with error";
  }
}
