/*!\file Tecplot_ExecutionMode.cc
  \brief Implementation of member functions of class Tecplot_Execution_Mode. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Tecplot_ExecutionMode.h"
#include "../MPI/MPI.h"

short Tecplot_Execution_Mode::OUTPUT_LEVEL = Detailed; // Set a detail level of output
short Tecplot_Execution_Mode::OUTPUT_ACCURACY = SinglePrec; // Set to write the output data in single precision

//! Set all flags to default values
// add all flag default values to this function
void Tecplot_Execution_Mode::SetDefaults(void){

  OUTPUT_LEVEL = Detailed; // Set a detail level of output
  OUTPUT_ACCURACY = SinglePrec; // Set to write the output data in single precision
}

//! Print the current execution mode
//  at the output stream
// \param [in] out_file the output stream
void Tecplot_Execution_Mode::Print_Info(std::ostream & out_file){

  // write what the output format is
  if (OUTPUT_LEVEL == Brief){
    out_file << "\n     -> Tecplot Output Level: " << "Brief";
  } else if (OUTPUT_LEVEL == Detailed){
    out_file << "\n     -> Tecplot Output Level: " << "Detailed";
  } else if (OUTPUT_LEVEL == Full){
    out_file << "\n     -> Tecplot Output Level: " << "Full";
  } // endif

  if (OUTPUT_ACCURACY == SinglePrec){
    out_file << "\n     -> Tecplot Precision: " << "Single";
  } else {
    out_file << "\n     -> Tecplot Precision: " << "Double";
  }
}

/*!
 * Broadcast the Tecplot_Execution_Mode variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Tecplot_Execution_Mode::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&OUTPUT_LEVEL,
 			1, 
 			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&OUTPUT_ACCURACY,
 			1, 
 			MPI::SHORT, 0);

#endif
}
