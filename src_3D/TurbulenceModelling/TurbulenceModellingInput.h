/* TurbulenceModellingInput.h:  Header file declaring turbulence model
                                input class. */

#ifndef _TURBULENCEMODEL_INPUT_INCLUDED
#define _TURBULENCEMODEL_INPUT_INCLUDED

// Include required C++ libraries.

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

// Include CFD and MPI header files.

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

#ifndef _MPI_INCLUDED
#include "../MPI/MPI.h"
#endif // _MPI_INCLUDED

#define TURBULENCEMODEL_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: Turbulence_Modelling_Input_Parameters
 *
 * @brief Input Parameters for the modelling of turbulence flows.
 *
 * This class defines and handles the input variables related to the
 * modelling and treatment of turbulent flows.
 *
 */
class Turbulence_Modelling_Input_Parameters{
 private:
 public:
  //@{ @name Indicator for near wall treatment for turbulent flows using FANS solver:
  int Near_Wall_Boundary_Treatment; 
  //! 0, 1,2 , automatic, wall function, low_Reynolds number.
  //@}

  //@{ @name Constructor and desctructor
  //! Constructor (assign default values).
  Turbulence_Modelling_Input_Parameters() {
    Near_Wall_Boundary_Treatment = 0; 
  }

  //! Destructor
  ~Turbulence_Modelling_Input_Parameters(void){ }
  //@}

  //@{ @name Other Member functions.
  //! Broadcast input parameters to all processors:
  void Broadcast(void);
  //! Parse next input line:
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);
  //! Check validity of specified input parameters:
  int Check_Inputs(void);
  //@}

  //@{ @name Input-output operators:
  friend ostream &operator << (ostream &out_file,
		               const Turbulence_Modelling_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Turbulence_Modelling_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  //@}

};

#endif // _TURBULENCEMODEL_INPUT_INCLUDED
