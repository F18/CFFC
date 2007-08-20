/* MPI.cc:  Subroutines and Variable for MPI (Message Passing Interface) Class. */

/* Include MPI header file. */

#ifndef _MPI_INCLUDED
#include "MPI.h"
#endif // _MPI_INCLUDED

/*********************************************************************
 * CFFC_MPI -- Create storage and initialize global MPI variables. *
 *********************************************************************/
int CFFC_MPI::Number_of_Processors = 0;
int CFFC_MPI::This_Processor_Number = 0;
