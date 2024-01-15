/*!\file CENO_ExecutionMode.cc
  \brief Initialize the flags that control the execution of CENO high-order reconstruction. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../../../src_2D/CFD/CFD.h"
#include "CENO_CFRS_ExecutionMode.h"

short CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE = OFF; // computation memory efficient mode
short CENO_CFRS_Execution_Mode::REDUCE_ORDER       = ON;  // recontruct to limited piecewise linear based on SI
//short CENO_CFRS_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING = ON;

//! Set all flags to default values
// add all flag default values to this function
void CENO_CFRS_Execution_Mode::SetDefaults(void){
  
  USE_PSEUDO_INVERSE = OFF; // computation memory efficient mode
  REDUCE_ORDER       = ON;  // recontruct to limited piecewise linear based on SI
  // CENO_APPLY_GEOMETRIC_WEIGHTING == ON
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void CENO_CFRS_Execution_Mode::Print_Info(std::ostream & out_file){

  // output execution mode
  if (USE_PSEUDO_INVERSE == OFF){
    out_file << "\n     -> Execution Mode: " << "Memory Efficient: Pseudo Inverse Not Used";
    out_file << "\n     -> Least-squares solver";
  } else {
    out_file << "\n     -> Execution Mode: " << "Speed Efficient: Pseudo Inverse Used";
  }

  if (REDUCE_ORDER == ON){
    out_file << "\n     -> Execution Mode: " << "Reduce to Limited Piecewise Linear reconstrcution based on SI.";
  } else {
    out_file << "\n     -> Execution Mode: " << "No reduction of order. Keep Unlimited High-Order reconstruction";
  }

//  // output geom weighting type
//  if (CENO_APPLY_GEOMETRIC_WEIGHTING == ON){
//    if (CENO_SQUARE_GEOM_WEIGHTING == ON){
//      out_file << "\n     -> Geom Weighting: " << "Inverse of squared distance";
//    } else {
//      out_file << "\n     -> Geom Weighting: " << "Inverse distance";
//    }
//  } else {
//    out_file << "\n     -> Geom Weighting: " << "Turned OFF";
//  }

}

/*!
 * Broadcast the CENO_CFRS_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
//void CENO_CFRS_Execution_Mode::Broadcast(void){
//#ifdef _MPI_VERSION
//  
//  MPI::COMM_WORLD.Bcast(&USE_CENO_ALGORITHM,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_PSEUDO_INVERSE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_DROP_ORDER,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_PADDING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_APPLY_GEOMETRIC_WEIGHTING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_SQUARE_GEOM_WEIGHTING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSIDER_WEIGHTS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_LAPACK_LEAST_SQUARES,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&HIGH_ORDER_MESSAGE_PASSING,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&CENO_VERBOSE,
// 			1, 
// 			MPI::SHORT, 0);
//  MPI::COMM_WORLD.Bcast(&Limiter,
// 			1, 
// 			MPI::SHORT, 0);  
//
//#endif
//}
