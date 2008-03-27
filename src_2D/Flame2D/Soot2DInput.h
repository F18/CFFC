/////////////////////////////////////////////////////////////////////
///
/// \file SootState.h
/// 
/// \author Marc R.J. Charest
/// 
/// \brief This header file contains the class definition for 
///        describing the soot inputs.
/// 
/////////////////////////////////////////////////////////////////////
#ifndef _SOOT2D_INPUT_INCLUDED 
#define _SOOT2D_INPUT_INCLUDED

// CFFC includes
#include "Soot2DState.h"

/////////////////////////////////////////////////////////////////////
/// Prototypes
/////////////////////////////////////////////////////////////////////
//! set parameter values
int setParamValue( int &param,             // parameter to set
		   const char *value,      // value to compare with list
		   const int command,      // the command value
		   const string *list,     // list of possible values
		   const int n );          // number of values

/////////////////////////////////////////////////////////////////////
/// CLASS DEFINITIONS
/////////////////////////////////////////////////////////////////////
/**
 * \class Soot2D_Input_Parameters
 *
 * Soot State Inputs class definition.
 *
 */
class Soot2D_Input_Parameters{
 public:

  //! flag whether mixture is sooting
  int sootModel;

  //! Default Constructor
  Soot2D_Input_Parameters() {
    sootModel = SOOT2D_NONE;
  };

  //! Accessors
  string sootModelName(void) const { return SOOT_MODEL_NAMES[sootModel]; };

  // output parameters
  void Output(ostream &out_file) const;

  // input file parser
  int Parse_Next_Input_Control_Parameter(char *code, char *value);

  // MPI communicators
  void Broadcast_Input_Parameters();

  // IO operators  
  friend ostream &operator << (ostream &out_file,
		               const Soot2D_Input_Parameters &IP);
  friend istream &operator >> (istream &in_file,
			       Soot2D_Input_Parameters &IP);
};

/////////////////////////////////////////////////////////////////////
/// INLINE FUNCTIONS
/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Soot2D_Input_Parameters :: Broadcast_Input_Parameters              
 *                                                                   
 * Broadcast the input parameters variables to all processors        
 * involved in the calculation from the primary processor using the  
 * MPI broadcast routine.                                            
 *********************************************************************/
inline void Soot2D_Input_Parameters :: Broadcast_Input_Parameters() {
#ifdef _MPI_VERSION
  MPI::COMM_WORLD.Bcast(&(sootModel), 1, MPI::INT, 0);
#endif
}

/*********************************************************************
 * Soot2D_Input_Parameters :: Parse_Next_Input_Control_Parameter     
 *                                                                   
 * Get the next input control parameter from the input file.         
 * Returns:                                                          
 *  - INVALID_INPUT_VALUE if code is valid but value is invalid      
 *  - INVALID_INPUT_CODE  if unknown code                            
 *********************************************************************/
inline int Soot2D_Input_Parameters :: 
Parse_Next_Input_Control_Parameter(char *code, char *value)
{
  int i_command = INVALID_INPUT_CODE;

  if (strcmp(code, "Soot_Model") == 0)  {
    i_command = setParamValue(sootModel, value, 9000, SOOT_MODEL_NAMES, NUM_SOOT_MODELS);

  } else {
    i_command = INVALID_INPUT_CODE;
  }

  return i_command;
  
}


/*********************************************************************
 * SNBCK_Input_Parameters :: Output                                  *
 *                                                                   *
 * Display Output Operator.                                          *
 *********************************************************************/
inline void Soot2D_Input_Parameters::Output(ostream& out) const{
  if (!sootModel) out << "\n  -> Soot Model   ====> OFF";
  else            out << "\n  -> Soot Model   ====> " << sootModelName();
}

/*********************************************************************
 * setParamValue                                                     *
 *                                                                   *
 * Set a parameter value.                                            *
 *********************************************************************/
inline int setParamValue( int &param,             // parameter to set
			  const char* value,      // value to compare with list
			  const int command,      // the command value
			  const string *list,     // list of possible values
			  const int n )  {        // number of values

  // search through the list
  for (int i=0; i<n; i++) {
    if (strcmp(value, list[i].c_str()) == 0) {
      param = i;
      return command;
    }
  }

  // not found, return invalid input
  return INVALID_INPUT_CODE;
}


#endif // _SOOT2D_INPUT_INCLUDED
