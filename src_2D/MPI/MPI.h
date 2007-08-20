/* MPI.h:  Header file defining required MPI (Message Passing Interface) 
           subroutines and macros. */

#ifndef _MPI_INCLUDED
#define _MPI_INCLUDED

/* Include CFFC CFD header file. */

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

/* Include MPI C++ Language Binding header file. */

#ifdef _MPI_VERSION
#ifdef _MPI_VERSION2
//#undef SEEK_SET
//#undef SEEK_CUR
//#undef SEEK_END
#include "mpi.h"
#else
#include "mpi++.h"
#endif
#endif

/* Define the classes. */

/********************************************************
 * Class: CFFC_MPI                                    *
 *                                                      *
 * Member functions                                     *
 *   Number_of_Processors  -- Return number of MPI      *
 *                            processors.               *
 *   This_Processor_Number -- Return MPI processor      *
 *                            number.                   *
 *                                                      *
 * Member operators                                     *
 *      M -- a MPI global variable                      *
 *                                                      *
 * cout << M; (output function)                         *
 * cin  >> M; (input function)                          *
 *                                                      *
 ********************************************************/
class CFFC_MPI{
  private:
  public:
    static int Number_of_Processors;   // Number of processors.
    static int This_Processor_Number;  // Processor number.
                                       // Made public so can access them.
		      
    /* Creation, copy, and assignment constructors. */
    //CFFC_MPI(void) { }
    // Use automatically generated constructor.

    /* Destructor. */
    // ~CFFC_MPI(void);
    // Use automatically generated destructor.

    /* Input-output operators. */

    friend ostream &operator << (ostream &out_file, const CFFC_MPI &mpi);
    friend istream &operator >> (istream &in_file, CFFC_MPI &mpi);

};

/*************************************************************
 * CFFC_MPI -- Input-output operators.                     *
 *************************************************************/
inline ostream &operator << (ostream &out_file, const CFFC_MPI &mpi) {
    return (out_file);
}

inline istream &operator >> (istream &in_file, CFFC_MPI &mpi) {
    return (in_file);
}

/********************************************************
 * MPI -- Inline functions.                             *
 ********************************************************/

// MPI Version.
inline string CFFC_Version_MPI() {
  int version, subversion; 
  char string1[128], string2[128];
#ifdef _MPI_VERSION
#if MPI2CPP_HAVE_MPI_GET_VERSION
  MPI::Get_version(version, subversion);
  strcpy(string1, "MPI Version, Built using Version ");
  sprintf(string2, "%d.%d", version, subversion);
  strcat(string2, " of Message-Passing Interface (MPI) Standard.");
  strcat(string1, string2);
  return (string(string1));
#else
  return ("MPI Version, Built using Version 1.xx of Message-Passing Interface (MPI) Standard.");
#endif
#endif
#ifdef _NO_MPI_VERSION
  return ("No Message-Passing Interface (MPI) Version.");
#endif
}

// Initialize MPI.
inline void CFFC_Initialize_MPI(int num_arg, char *arg_ptr[]) {
#ifdef _MPI_VERSION
  MPI::Init(num_arg, arg_ptr);
  CFFC_MPI::Number_of_Processors  = MPI::COMM_WORLD.Get_size();
  CFFC_MPI::This_Processor_Number = MPI::COMM_WORLD.Get_rank();
#else
  CFFC_MPI::Number_of_Processors  = 1;
  CFFC_MPI::This_Processor_Number = 0;
#endif
}

// Finalize MPI.
inline void CFFC_Finalize_MPI() {
#ifdef _MPI_VERSION
  MPI::Finalize();
#endif
}

// Primary MPI Processor.
inline int CFFC_Primary_MPI_Processor() {
  return(CFFC_MPI::This_Processor_Number == 0);
}

// MPI Barrier for all processors.
inline void CFFC_Barrier_MPI() {
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Barrier();
#endif
}

// Broadcast integers to all processors.
inline void CFFC_Broadcast_MPI(int *buffer, const int buffer_size) {
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::INT, 0);
#endif
}

// Broadcast doubles to all processors.
inline void CFFC_Broadcast_MPI(double *buffer, const int buffer_size) {
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);
#endif
}

// Returns AND of all integers sent by each processor.
inline int CFFC_AND_MPI(const int send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::LAND);
#endif
  return (receive_value);
}

// Returns OR of all integers sent by each processor.
inline int CFFC_OR_MPI(const int send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::LOR);
#endif
  return (receive_value);
}

// Returns minimum value of all integers sent by each processor.
inline int CFFC_Minimum_MPI(const int send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::MIN);
#endif
  return (receive_value);
}

// Returns minimum value of all doubles sent by each processor.
inline double CFFC_Minimum_MPI(double &send_value) {
  double receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::DOUBLE, MPI::MIN);
#endif
  return (receive_value);
}

// Returns maximum value of all integers sent by each processor.
inline int CFFC_Maximum_MPI(const int &send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::MAX);
#endif
  return (receive_value);
}

// Returns maximum value of all doubles sent by each processor.
inline double CFFC_Maximum_MPI(double &send_value) {
  double receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::DOUBLE, MPI::MAX);
#endif
  return (receive_value);
}

// Returns summation of all integers sent by each processor.
inline int CFFC_Summation_MPI(const int send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::SUM);
#endif
  return (receive_value);
}

// Returns summation of all doubles sent by each processor.
inline double CFFC_Summation_MPI(double &send_value) {
  double receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::DOUBLE, MPI::SUM);
#endif
  return (receive_value);
}

// Returns summation of all state values sent by each processor.
template <class T> inline T CFFC_Summation_MPI(T &send_value) {
  T receive_value(send_value);
#ifdef _MPI_VERSION
  for (int n = 1; n < send_value.NumVar(); n++)
  MPI::COMM_WORLD.Allreduce(&send_value[n], &receive_value[n], 1, MPI::DOUBLE, MPI::SUM);
#endif
  return (receive_value);
}

// Returns product of all integers sent by each processor.
inline int CFFC_Product_MPI(const int send_value) {
  int receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::INT, MPI::PROD);
#endif
  return (receive_value);
}

// Returns product of all doubles sent by each processor.
inline double CFFC_Product_MPI(double &send_value) {
  double receive_value;
  receive_value = send_value;
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Allreduce(&send_value, &receive_value, 1, MPI::DOUBLE, MPI::PROD);
#endif
  return (receive_value);
}

/********************************************************
 * MPI -- Define MPI structures and classes.            *
 ********************************************************/

/********************************************************
 * MPI -- External subroutines.                         *
 ********************************************************/

#endif /* _MPI_INCLUDED  */
