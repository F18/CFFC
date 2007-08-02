/* MPI.cc:  Subroutines and Variable for MPI (Message Passing Interface) Class. */

/* Include MPI header file. */

#ifndef _MPI_INCLUDED
#include "MPI.h"
#endif // _MPI_INCLUDED

/*********************************************************************
 * CFDkit_MPI -- Create storage and initialize global MPI variables. *
 *********************************************************************/
int CFDkit_MPI::Number_of_Processors = 0;
int CFDkit_MPI::This_Processor_Number = 0;
