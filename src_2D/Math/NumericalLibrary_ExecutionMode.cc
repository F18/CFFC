/*!\file NumericalLibrary_ExecutionMode.cc
  \brief Source file defining the values of the numerical parameters/flags defined in NumericalLibrary_ExecutionMode.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "NumericalLibrary_ExecutionMode.h"


// NumericalLibrary_Execution_Mode class

/*!
 * This value is used to set the number of function evaluations allowed
 * during the adaptive integration.
 * It's the stopping criterion for not converged integrals.
 */
long int NumericalLibrary_Execution_Mode::Max_Function_Evaluations = 8000000;

/*!
 * This value is used to set the number of samples for the Monte Carlo integration.
 */
long int NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples = 200000;

/*!
 * This value is used to set the minimum number of refinement levels in 
 * the adaptive integration algorithm.
 */
short NumericalLibrary_Execution_Mode::Adaptive_Integration_Minimum_Refinement_Levels = 0;

/*!
 * This value controls the error output.
 */
bool NumericalLibrary_Execution_Mode::Output_Error_Messages = true;

// Set default epsilon values
long int NumericalLibrary_Execution_Mode::Max_Function_Evaluations_default = 
  NumericalLibrary_Execution_Mode::Max_Function_Evaluations;
long int NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples_default = 
  NumericalLibrary_Execution_Mode::Number_Monte_Carlo_Samples;
short NumericalLibrary_Execution_Mode::Adaptive_Integration_Minimum_Refinement_Levels_default = 
  NumericalLibrary_Execution_Mode::Adaptive_Integration_Minimum_Refinement_Levels;

/*! Print the current execution mode
 *  at the output stream
 * \param [in] out_file the output stream
 */
void NumericalLibrary_Execution_Mode::Print_Info(std::ostream & out_file){

  // output Max_Function_Evaluations
  if (Max_Function_Evaluations != Max_Function_Evaluations_default){
    out_file << "\n     -> Max allowed of function evaluations during integration: " << Max_Function_Evaluations
	     << "  (default value = " << Max_Function_Evaluations_default << ")";
  }

  // output Number_Monte_Carlo_Samples
  if (Number_Monte_Carlo_Samples != Number_Monte_Carlo_Samples_default){
    out_file << "\n     -> Number of Monte Carlo samples during integration: " << Number_Monte_Carlo_Samples
	     << "  (default value = " << Number_Monte_Carlo_Samples_default << ")";
  }

  // output Adaptive_Integration_Minimum_Refinement_Levels
  if (Adaptive_Integration_Minimum_Refinement_Levels != Adaptive_Integration_Minimum_Refinement_Levels_default){
    out_file << "\n     -> Minimum number of adaptive integration levels: " 
	     << Adaptive_Integration_Minimum_Refinement_Levels
	     << "  (default value = " << Adaptive_Integration_Minimum_Refinement_Levels_default << ")";
  }
  
}

/*!
 * Set default parameter values
 */
void NumericalLibrary_Execution_Mode::SetDefaults(void){

  Max_Function_Evaluations = Max_Function_Evaluations_default;
  Number_Monte_Carlo_Samples = Number_Monte_Carlo_Samples_default;
  Adaptive_Integration_Minimum_Refinement_Levels = Adaptive_Integration_Minimum_Refinement_Levels_default;

}

/*!
 * Broadcast the NumericalLibrary_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 *
 */
void NumericalLibrary_Execution_Mode::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&Max_Function_Evaluations,
			1, 
			MPI::LONG, 0);
  MPI::COMM_WORLD.Bcast(&Number_Monte_Carlo_Samples,
			1, 
			MPI::LONG, 0);
  MPI::COMM_WORLD.Bcast(&Adaptive_Integration_Minimum_Refinement_Levels,
			1, 
			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&Output_Error_Messages,
			1, 
			MPI::BOOL, 0);
#endif
}
