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
vector<int> HighOrder2D_Input::OrdersOfReconstruction(1,-1); //!< Set to no high-order reconstruction.
int HighOrder2D_Input::NumberOfAuxiliaryReconstructions(0); //!< Set to ZERO. No auxiliary reconstructions.
int HighOrder2D_Input::NumberOfHighOrderReconstructions(0); //!< Set to ZERO. No high-order reconstruction.
int HighOrder2D_Input::MaxReconstructionOrder(0);		  //!< Set to ZERO (i.e piecewise constant ).


// ===  Member functions ===
/*!
 * Reset class variables to default values.
 */
void HighOrder2D_Input::SetDefaults(void){
  OrdersOfReconstruction = vector<int>(1,-1); // Set to no high-order reconstruction.
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
 */
void HighOrder2D_Input::Broadcast(void){
#ifdef _MPI_VERSION
  
  int i;
  int buffer_size;
  int *buffer;

  // On primary CPU
  if (CFFC_Primary_MPI_Processor()) {
    // calculate buffer_size
    buffer_size = 3 + OrdersOfReconstruction.size();

    // allocate and load buffer
    buffer = new int [buffer_size];

    buffer[0] = NumberOfHighOrderReconstructions;
    buffer[1] = NumberOfAuxiliaryReconstructions;
    buffer[2] = MaxReconstructionOrder;

    for (i = 3; i<buffer_size; ++i){
      buffer[i] = OrdersOfReconstruction[i-3];
    }
  }

  // Broadcast the buffer size
  MPI::COMM_WORLD.Bcast(&buffer_size, 1, MPI::INT, 0);
 
  /* On non-primary MPI processors, allocate 
     memory for the buffer. */
  if (!CFFC_Primary_MPI_Processor()) {
    buffer = new int [buffer_size];
  } /* endif */

  // Broadcast the buffer
  MPI::COMM_WORLD.Bcast(&buffer[0], buffer_size, MPI::INT, 0);

  // On non-primary MPI processors, unpack the buffer.
  if (!CFFC_Primary_MPI_Processor()) {
    NumberOfHighOrderReconstructions = buffer[0];
    NumberOfAuxiliaryReconstructions = buffer[1];   
    MaxReconstructionOrder = buffer[2];

    OrdersOfReconstruction.assign(NumberOfHighOrderReconstructions,-1);
    for (i = 3; i<buffer_size; ++i){
      OrdersOfReconstruction[i-3] = buffer[i];
    }
  }
  
#endif
}
