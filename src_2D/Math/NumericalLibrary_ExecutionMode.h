/*!\file NumericalLibrary_ExecutionMode.h
  \brief Definition of flags that control the execution of routines in the Numerical Library. */

#ifndef _NUMERICALLIBRARY_EXECUTIONMODE_INCLUDED
#define _NUMERICALLIBRARY_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"


/*!
 * \class NumericalLibrary_Execution_Mode
 *
 * @brief Definition of flags that control the execution of routines in the Numerical Library
 *
 */
class NumericalLibrary_Execution_Mode{
  
public:

  /*! This parameter controls the maximum number of function evaluations in the integration routines.\n
      -----------------------------------------------------------------------------------------------  */
  static long int Max_Function_Evaluations;
  
  /*! This parameter controls the number of samples in the Monte Carlo integration routines.\n
      -------------------------------------------------------------------------------------  */
  static long int Number_Monte_Carlo_Samples;

  /*! This parameter controls the minimum number of refinement levels in the adaptive integration routines.\n
      ---------------------------------------------------------------------------------------------------  */
  static short Adaptive_Integration_Minimum_Refinement_Levels;

  /*! This parameter controls whether the error messages are output. \n
   *  DON'T TURN OFF THE ERROR OUTPUT UNLESS YOU HAVE STRONG REASONS!!!!
      ---------------------------------------------------------------------------------------------------  */
  static short Output_Error_Messages;
  
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  // Set all parameters to default values
  static void SetDefaults(void);

  static void Broadcast(void);

protected:
  NumericalLibrary_Execution_Mode(void);   //!< Private default constructor
  NumericalLibrary_Execution_Mode(const NumericalLibrary_Execution_Mode&); //!< Private copy constructor
  NumericalLibrary_Execution_Mode& operator=(const NumericalLibrary_Execution_Mode&); //!< Private assignment operator

  //! @name Copy of class parameters that cannot be modified at runtime
  //@{   
  static long int Max_Function_Evaluations_default;
  static long int Number_Monte_Carlo_Samples_default;
  static short Adaptive_Integration_Minimum_Refinement_Levels_default;
  //@}

};

//! Parse the input control parameters for 
//  settings related to NumericalLibrary_Execution_Mode class
template<class Input_Parameters_Type> inline
void NumericalLibrary_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP,
									 int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Integration_Maximum_Function_Evaluation") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Max_Function_Evaluations;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Number_Monte_Carlo_Samples") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Number_Monte_Carlo_Samples;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Adaptive_Integration_Minimum_Refinement_Levels") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Adaptive_Integration_Minimum_Refinement_Levels;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Numerical_Library_Output_Error_Messages") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 || strcmp(IP.Next_Control_Parameter, "YES") == 0 ){
      Output_Error_Messages = ON;
      i_command = 0;
    } else if ( strcmp(IP.Next_Control_Parameter, "No") == 0 || strcmp(IP.Next_Control_Parameter, "NO") == 0 ){
      Output_Error_Messages = OFF;
      i_command = 0;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif
