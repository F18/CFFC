/*!\file CENO_ExecutionMode.cc
  \brief Initialize the flags that control the execution of CENO high-order reconstruction. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../CFD/CFD.h"
#include "CENO_ExecutionMode.h"
#include "../MPI/MPI.h"

short CENO_Execution_Mode::USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
short CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
short CENO_Execution_Mode::CENO_DROP_ORDER = ON; // produce monotone solutions
short CENO_Execution_Mode::CENO_PADDING = OFF; // no padding
short CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING = ON; // apply geometric weighting
short CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
short CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights
short CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface
short CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF; // use enough layers of ghost cells
short CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF; // compute SI with all neighbours
short CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = ON; //use additional constraints
short CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = OFF; // don't extend the stencil
short CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES = ON; // use Lapack least-squares subroutine
int CENO_Execution_Mode::Limiter = LIMITER_VANLEER;
short CENO_Execution_Mode::IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = OFF; // don't ignore curved boundaries
short CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = ON; // use the smoothness indicator for CENO AMR

//! Set all flags to default values
// add all flag default values to this function
void CENO_Execution_Mode::SetDefaults(void){
  
  USE_CENO_ALGORITHM = OFF; // CENO scheme is not used
  CENO_SPEED_EFFICIENT = ON; // computation time efficient mode
  CENO_DROP_ORDER = ON; // produce monotone solutions
  CENO_PADDING = OFF; // no padding
  CENO_APPLY_GEOMETRIC_WEIGHTING = ON; // apply geometric weighting
  CENO_SQUARE_GEOM_WEIGHTING = OFF; // use 1.0/fabs(Distance)
  CENO_CONSIDER_WEIGHTS = OFF;	// computation of smoothness indicator without weights
  FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON; // try to use the PWC at interface
  CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF; // use enough layers of ghost cells to compute everything that is needed
  CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF; // compute SI with all neighbours
  CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = ON; // use additional constraints
  CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = OFF; // don't extend the stencil
  USE_LAPACK_LEAST_SQUARES = ON; // use Lapack least-squares subroutine
  IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = OFF; // don't ignore curved boundaries
  USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = ON; // use the smoothness indicator for CENO AMR
  Limiter = LIMITER_VANLEER;
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void CENO_Execution_Mode::Print_Info(std::ostream & out_file){

  // output execution mode
  if (CENO_SPEED_EFFICIENT == ON){
    out_file << "\n     -> Execution Mode: " << "Speed Efficient";
  } else {
    out_file << "\n     -> Execution Mode: " << "Memory Efficient";
  }

  // output subroutine used to solver the least-squares problem
  if (CENO_SPEED_EFFICIENT == OFF){
    if (USE_LAPACK_LEAST_SQUARES){
      out_file << "\n     -> Least-squares solver: " << "Lapack routine";
    } else {
      out_file << "\n     -> Least-squares solver: " << "Internal routine";
    }
  }

  // output monotonicity mode
  if (CENO_DROP_ORDER == ON){
    out_file << "\n     -> Monotonicity Mode: " << "Yes (drop order)";
  } else {
    out_file << "\n     -> Monotonicity Mode: " << "No (don't drop order)";
  }

  // output padding mode
  if (CENO_PADDING == ON){
    out_file << "\n     -> Cell Padding: " << "Yes (flag adjacent cells too)";
  } else {
    out_file << "\n     -> Cell Padding: " << "No";
  }

  // output geom weighting type
  if (CENO_APPLY_GEOMETRIC_WEIGHTING == ON){
    if (CENO_SQUARE_GEOM_WEIGHTING == ON){
      out_file << "\n     -> Geom Weighting: " << "Inverse of squared distance";
    } else {
      out_file << "\n     -> Geom Weighting: " << "Inverse distance";
    }
  } else {
    out_file << "\n     -> Geom Weighting: " << "Turned OFF";
  }

  // output interface behavior
  if (FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE == ON){
    out_file << "\n     -> Negative interface solutions: " << "Force with PWC";
  } else {
    out_file << "\n     -> Negative interface solutions: " << "Exit with error";
  }

  // output message passing if ON
  if (CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING == ON){
    out_file << "\n     -> CENO for ghost cells: " << "With message passing";
  }

  // output how many neighbours are used for calculation of smoothness indicator
  if (CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS == ON){
    out_file << "\n     -> CENO Smoothness Indicator: " << "With first neighbours only";
  } else {
    out_file << "\n     -> CENO Smoothness Indicator: " << "With all neighbours";
  }

  // output behaviour in for cells with constrained reconstruction
  if (CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS == ON && 
      CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL == OFF ){
    out_file << "\n     -> CENO Additional Constraints: " << "Yes";
    out_file << "\n     -> CENO Extended Biased Stencil: " << "No";
  } else {
    out_file << "\n     -> CENO Additional Constraints: " << "No";
    out_file << "\n     -> CENO Extended Biased Stencil: " << "Yes";
  }

  // output accuracy assessment behaviour close to curved boundaries
  if ( IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT == ON ){
    out_file << "\n     -> Accuracy Measurement Near Curved Boundaries: " << "Ignore cells";
  } else {		   
    out_file << "\n     -> Accuracy Measurement Near Curved Boundaries: " << "Include cells";
  }

  // output high-order AMR criteria
  if ( USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA == ON ){
    out_file << "\n     -> High-order AMR: " << "Use Smoothness Indicator";
  } else {
    out_file << "\n     -> High-order AMR: " << "DON'T Use Smoothness Indicator";
  }
}

/*!
 * Broadcast the CENO_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void CENO_Execution_Mode::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&USE_CENO_ALGORITHM,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_SPEED_EFFICIENT,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_DROP_ORDER,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_PADDING,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_APPLY_GEOMETRIC_WEIGHTING,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_SQUARE_GEOM_WEIGHTING,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_CONSIDER_WEIGHTS,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&USE_LAPACK_LEAST_SQUARES,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&Limiter,
 			1, 
 			MPI::SHORT, 0);  
#endif
}
