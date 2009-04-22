/*!\file AccuracyAssessment_ExecutionMode.h
  \brief Definition of flags that control the execution of the accuracy assessment framework. */

#ifndef _ACCURACY_ASSESSMENT_EXECUTIONMODE_INCLUDED
#define _ACCURACY_ASSESSMENT_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

/*!
 * Class: AccuracyAssessment_Execution_Mode
 *
 * @brief Definition of flags that control the execution of the accuracy assessment framework
 *
 */
class AccuracyAssessment_Execution_Mode{
  
public:
  //! Data type defining different ways of assessing the accuracy.
  enum Accuracy_Assessment_Type { Based_On_Exact_Solution,    //!< Assess the accuracy using an exact solution
				  Based_On_Entropy_Variation, /*!< Assess the accuracy based on the 
								production of entropy. Useful for inviscid flows. */
				  Based_On_Lift_And_Drag_Coefficients  /*!< Assess the accuracy based on the
									 values of the lift and drag coefficients. */
  };

  static void SetDefaults(void);

  //@{ @name Field access.
  static const short& Method(void) {return Accuracy_Assessment_Method; }
  static const unsigned int & Exact_Digits(void) {return Accuracy_Assessment_Exact_Digits; }
  static const unsigned int & Assessment_Parameter(void) {return Accuracy_Assessment_Parameter; }
  static const unsigned int & Assessment_Frequency(void) {return Accuracy_Assessment_Frequency; }
  //@}

  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  static void Broadcast(void);

protected:
  AccuracyAssessment_Execution_Mode(void);   //!< Private default constructor
  AccuracyAssessment_Execution_Mode(const AccuracyAssessment_Execution_Mode&); //!< Private copy constructor
  AccuracyAssessment_Execution_Mode& operator=(const AccuracyAssessment_Execution_Mode&); //!< Private assignment operator

  static short Accuracy_Assessment_Method;
  static unsigned int Accuracy_Assessment_Exact_Digits;
  static unsigned int Accuracy_Assessment_Parameter;
  static unsigned int Accuracy_Assessment_Frequency; //!< The frequency for assessing errors in time-dependent problems
};

//! Parse the input control parameters for 
//  settings related to AccuracyAssessment_Execution_Mode class
template<class Input_Parameters_Type> inline
void AccuracyAssessment_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Method") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Exact_Solution") == 0 ){
      Accuracy_Assessment_Method = Based_On_Exact_Solution;
    } else if ( strcmp(IP.Next_Control_Parameter, "Entropy_Variation") == 0 ) {
      Accuracy_Assessment_Method = Based_On_Entropy_Variation;
    } else if ( strcmp(IP.Next_Control_Parameter, "Calculate_Lift_And_Drag") == 0 ) {
      Accuracy_Assessment_Method = Based_On_Lift_And_Drag_Coefficients;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Exact_Digits") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Accuracy_Assessment_Exact_Digits;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (Accuracy_Assessment_Exact_Digits < 0) i_command = INVALID_INPUT_VALUE;
      
  } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Parameter") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Accuracy_Assessment_Parameter;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Frequency") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Accuracy_Assessment_Frequency;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (Accuracy_Assessment_Frequency < 0) i_command = INVALID_INPUT_VALUE;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif	// _ACCURACY_ASSESSMENT_EXECUTIONMODE_INCLUDED
