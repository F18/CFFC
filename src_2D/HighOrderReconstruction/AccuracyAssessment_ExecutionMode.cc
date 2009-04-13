/*!\file AccuracyAssessment_ExecutionMode.cc
  \brief Implementation of member functions of class AccuracyAssessment_Execution_Mode. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AccuracyAssessment_ExecutionMode.h"
#include "../MPI/MPI.h"

short AccuracyAssessment_Execution_Mode::Accuracy_Assessment_Method = Based_On_Exact_Solution; /* Set assessment of accuracy
												  with the exact solution */
unsigned int AccuracyAssessment_Execution_Mode::Accuracy_Assessment_Exact_Digits = 10;
unsigned int AccuracyAssessment_Execution_Mode::Accuracy_Assessment_Parameter = 1; // Use the first parameter (e.g. density)


//! Set all flags to default values
// add all flag default values to this function
void AccuracyAssessment_Execution_Mode::SetDefaults(void){

  Accuracy_Assessment_Method = Based_On_Exact_Solution; /* Set assessment of accuracy
							   with the exact solution */
  Accuracy_Assessment_Exact_Digits = 10;
  Accuracy_Assessment_Parameter = 1; // Use the first parameter (e.g. density)

}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void AccuracyAssessment_Execution_Mode::Print_Info(std::ostream & out_file){

  // output execution mode
  switch (Accuracy_Assessment_Method){
  case Based_On_Exact_Solution:
    out_file << "\n     -> Accuracy Assessment Method: " << "Relative to Exact Solution";
    break;
  case Based_On_Entropy_Variation:
    out_file << "\n     -> Accuracy Assessment Method: " << "Based on Entropy Variation";
    break;
  case Based_On_Lift_And_Drag_Coefficients:
    out_file << "\n     -> Accuracy Assessment Method: " << "Calculate Lift and Drag Coefficients";
    break;
  }

  // output number of exact digits
  out_file << "\n     -> Accuracy Assessment Exact Digits: " << Accuracy_Assessment_Exact_Digits;

  // output parameter
  if ( Accuracy_Assessment_Method == Based_On_Exact_Solution){
    out_file << "\n     -> Accuracy Assessment Parameter: " << Accuracy_Assessment_Parameter;
  }
}

/*!
 * Broadcast the AccuracyAssessment_Execution_Mode variables
 * to all processors associated with the specified communicator
 * from the specified processor using the MPI broadcast routine.
 */
void AccuracyAssessment_Execution_Mode::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&Accuracy_Assessment_Method,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&Accuracy_Assessment_Parameter,
 			1, 
 			MPI::INT, 0);
  MPI::COMM_WORLD.Bcast(&Accuracy_Assessment_Exact_Digits,
 			1, 
 			MPI::INT, 0);

#endif
}
