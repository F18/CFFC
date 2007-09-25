/* AMRInput.h:  Header file declaring AMR input class. */

#ifndef _AMRINPUT_INCLUDED
#define _AMRINPUT_INCLUDED

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

#define AMR_INPUT_PARAMETER_LENGTH 256

/*!
 * Class: AMR_Input_Parameters
 *
 * @brief Input Parameters for AMR Algorithm.
 *
 * This class defines and handles the input variables related to the
 * parallel implementation, domain decomposition, and block-based
 * adaptive mesh refinement (AMR).
 *
 */
class AMR_Input_Parameters{
 private:
 public:
  //@{ @name Multi-block solution-adaption and parallel domain decomposition input parameters:
  int Number_of_Processors, 
      Number_of_Blocks_Per_Processor;
  //@}

  //@{ @name AMR input parameters:
  //! Dynamic AMR flag:
  int AMR;
  //! Dynamic AMR frequency:
  int AMR_Frequency;
  //! Number of initial mesh refinements:
  int Number_of_Initial_Mesh_Refinements;
  //! Number of uniform mesh refinements:
  int Number_of_Uniform_Mesh_Refinements;
  //! Number of boundary mesh refinements:
  int Number_of_Boundary_Mesh_Refinements;
  //! Maximum number of refinement levels:
  int Maximum_Refinement_Level;
  //! Minimum number of refinement levels:
  int Minimum_Refinement_Level;
  //! Threshold for refinement of mesh:
  double Threshold_for_Refinement;
  //! Threshold for coarsening of mesh:
  double Threshold_for_Coarsening;
  //! Number of refinement criteria:
  int Number_of_Refinement_Criteria;
  //@}

  //@{ @name Morton ordering input parameters:
  int Morton;
  int Morton_Reordering_Frequency;
  //@}

  //@{ @name Constructor and desctructor
  //! Constructor (assign default values).
  AMR_Input_Parameters() {
     // Multi-block solution-adaption and parallel domain decomposition input parameters:
     Number_of_Processors = CFFC_MPI::Number_of_Processors;
     Number_of_Blocks_Per_Processor = 100;  
     // AMR input parameters:
     AMR = 0;
     AMR_Frequency = 100;
     Number_of_Initial_Mesh_Refinements = 0;
     Number_of_Uniform_Mesh_Refinements = 0;
     Number_of_Boundary_Mesh_Refinements = 0;
     Maximum_Refinement_Level = 100;
     Minimum_Refinement_Level = 1;
     Threshold_for_Refinement = 0.50;
     Threshold_for_Coarsening = 0.10;
     Number_of_Refinement_Criteria = 1;
     // Morton Ordering parameters:
     Morton = 0;
     Morton_Reordering_Frequency = 1000;
  }

  //! Destructor
  ~AMR_Input_Parameters(void) {}
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
		               const AMR_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       AMR_Input_Parameters &IP);
  void Output(ostream &out_file) const;
  //@}

};

#endif // _AMRINPUT_INCLUDED
