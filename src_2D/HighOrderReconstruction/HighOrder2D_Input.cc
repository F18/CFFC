/*!\file HighOrder2D_Input.cc
  \brief Source file initializing/implementing member variables/functions of class HighOrder2D_Input. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "HighOrder2D_Input.h"
#include "../Utilities/Utilities.h"

// ===  Static member variables ===
vector<short int> HighOrder2D_Input::OrdersOfReconstruction(1,-1); //!< Set to no high-order reconstruction.
short int HighOrder2D_Input::NumberOfAuxiliaryReconstructions(0); //!< Set to ZERO. No auxiliary reconstructions.
short int HighOrder2D_Input::NumberOfHighOrderReconstructions(0); //!< Set to ZERO. No high-order reconstruction.
short int HighOrder2D_Input::MaxReconstructionOrder(0);		  //!< Set to ZERO (i.e piecewise constant ).


// ===  Member functions ===
/*!
 * Reset class variables to default values.
 */
void HighOrder2D_Input::SetDefaults(void){
  OrdersOfReconstruction = vector<short int>(1,-1); // Set to no high-order reconstruction.
  NumberOfAuxiliaryReconstructions = 0; // Set to ZERO. No auxiliary reconstructions.
  NumberOfHighOrderReconstructions = 0; // Set to ZERO. No high-order reconstruction.
  MaxReconstructionOrder = 0;		// Set to ZERO (i.e piecewise constant ).
}

/*!
 * Print the relevant parameters of the HighOrder2D_Input
 * class to the provided output stream.
 */
void HighOrder2D_Input::Print_Info(std::ostream & out_file){

  if (NumberOfAuxiliaryReconstructions > 0){
    out_file << "\n  -> Auxiliary reconstructions: " << NumberOfAuxiliaryReconstructions;
    for (int i=1; i< OrdersOfReconstruction.size(); ++i){
      out_file << "\n     -> Reconstruction Order: " << OrdersOfReconstruction[i];
    }
  }

}

/*!
 * Broadcast the HighOrder2D_Input variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 *
 * \todo Check this subroutine!
 */
void HighOrder2D_Input::Broadcast(void){
#ifdef _MPI_VERSION
  
  int i;
  short int buffer_size;
  short int *buffer;

  // On primary CPU
  if (CFFC_Primary_MPI_Processor()) {
    // calculate buffer_size
    buffer_size = 2 + OrdersOfReconstruction.size();

    // allocate and load buffer
    buffer = new short int [buffer_size];

    buffer[0] = NumberOfHighOrderReconstructions;
    buffer[1] = NumberOfAuxiliaryReconstructions;

    for (i = 2; i<buffer_size; ++i){
      buffer[i] = OrdersOfReconstruction[i-2];
    }
  }

  // Broadcast the buffer size
  MPI::COMM_WORLD.Bcast(&buffer_size, 1, MPI::SHORT, 0);
  
  /* On non-primary MPI processors, allocate 
     memory for the buffer. */
  if (!CFFC_Primary_MPI_Processor()) {
    buffer = new short int [buffer_size];
  } /* endif */


  // Broadcast the buffer
  MPI::COMM_WORLD.Bcast(&buffer, buffer_size, MPI::SHORT, 0);


  // On non-primary MPI processors, unpack the buffer.
  if (!CFFC_Primary_MPI_Processor()) {
    NumberOfHighOrderReconstructions = buffer[0];
    NumberOfAuxiliaryReconstructions = buffer[1];    

    OrdersOfReconstruction.assign(NumberOfHighOrderReconstructions,-1);
    for (i = 2; i<buffer_size; ++i){
      OrdersOfReconstruction[i-2] = buffer[i];
    }
  }
  
#endif
}
