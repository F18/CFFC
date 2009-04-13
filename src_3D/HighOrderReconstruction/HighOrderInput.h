/* HighOrderInput.h:  Header file defining the 3D input data class. */

#ifndef _HIGHORDER_INPUT_INCLUDED
#define _HIGHORDER_INPUT_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include required CFFC header files. */

/* Include the CENO reconstruction execution mode header file. */
#ifndef _CENO_EXECUTION_MODE_INCLUDED
#include "../HighOrderReconstruction/CENO_ExecutionMode.h"
#endif // _GRID3D_HO_EXECUTIONMODE_INCLUDED

#ifndef _CENO_TOLERANCES_INCLUDED
#include "../HighOrderReconstruction/CENO_Tolerances.h"
#endif // _CENO_TOLERANCES

#define HIGHORDER_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: HighOrder_Input_Parameters
 *
 * @brief Input Parameters for High Order CENO Reconstruction.
 *
 * This class defines and handles the input variables related to the
 * the 3D high order CENO reconstruction.
 *
 */
class HighOrder_Input_Parameters{
  public:

  //@{ @name Constructors and desctructors:
  //! Constructor
  HighOrder_Input_Parameters(void){}

  //! Destructor
  ~HighOrder_Input_Parameters(void){}
  //@}

  //! Parse next input line
  int Parse_Next_Input_Control_Parameter(char *code, stringstream &value);

//  //@{ @name Input-output operators:
//  friend ostream &operator << (ostream &out_file,
//			       const HighOrder_Input_Parameters &IP);
//  friend istream &operator >> (istream &in_file,
//			       HighOrder_Input_Parameters &IP);
//  void Output(ostream &out_file) const;
//  //@}

};

#endif // _HIGHORDER_INPUT_INCLUDED
