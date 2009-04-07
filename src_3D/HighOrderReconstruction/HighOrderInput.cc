/* HighOrderInput.cc: Definition of HighOrder_Input_Parameters class member functions. */

/* Include the HighOrderInput header file. */

#ifndef _HIGHORDER_INPUT_INCLUDED
#include "HighOrderInput.h"
#endif // _HIGHORDER_INPUT_INCLUDED

/* Define member functions. */

/*********************************************************************************
 * HighOrder_Input_Parameters::Parse_Next_Input_Control_Parameter - Parse input.    *
 *********************************************************************************/
int HighOrder_Input_Parameters::Parse_Next_Input_Control_Parameter(char *code, 
								   stringstream &value) {

// Returns:
//  - INVALID_INPUT_VALUE if code is valid but value is invalid
//  - INVALID_INPUT_CODE  if unknown code

  int i_command = INVALID_INPUT_CODE;
  string value_string;


  /* ----- Parse input to set the static variables in the CENO_Execution_Mode class ----- */

  if (strcmp(code, "CENO_Execution_Mode") == 0) {
    value >> value_string;
    if (value_string == "Speed_Efficient"){
      CENO_Execution_Mode::CENO_SPEED_EFFICIENT = ON;
    } else if (value_string == "Memory_Efficient") {
      CENO_Execution_Mode::CENO_SPEED_EFFICIENT = OFF;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return i_command;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Drop_Order") == 0) {
    value >> value_string;
    if ( value_string == "Yes") {
      CENO_Execution_Mode::CENO_DROP_ORDER = ON;      
    } else {
      CENO_Execution_Mode::CENO_DROP_ORDER = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Padding") == 0) {
    value >> value_string;
    if ( strcmp(code, "Yes") == 0 ){
      CENO_Execution_Mode::CENO_PADDING = ON;
    } else {
      CENO_Execution_Mode::CENO_PADDING = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Apply_Geom_Weighting") == 0) {
    value >> value_string;
    if ( value_string == "Yes") {
      CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING = ON;
    } else {
      CENO_Execution_Mode::CENO_APPLY_GEOMETRIC_WEIGHTING = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Use_Geom_Weighting_For_Smoothness_Analysis") == 0) {
    value >> value_string;
    if (value_string == "Yes") {
      CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS = ON;
    } else {
      CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Square_Geom_Weighting") == 0) {
    value >> value_string;
    if (value_string == "Yes") {
      CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING = ON;
    } else {
      CENO_Execution_Mode::CENO_SQUARE_GEOM_WEIGHTING = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_Force_With_PWC") == 0) {
    value >> value_string;
    if (value_string == "Yes") {    
      CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON;
    } else {
      CENO_Execution_Mode::FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = OFF;
    }
    i_command = 0;
  } else if (strcmp(code, "CENO_With_Message_Passing") == 0) {
    value >> value_string;
    if (value_string == "Yes") {
      CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = ON;
      throw runtime_error("CENO_Execution_Mode ERROR! CENO reconstruction with message passing has not been implemented yet!");
    } else {
      CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "CENO_Smoothness_Indicator") == 0) {
    value >> value_string;
    if (value_string == "Use_First_Neighbours") {
      CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;
    } else if (value_string == "Use_All_Neighbours") {
      CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "Least_Squares_Solver") == 0) {
    value >> value_string;
    if ( value_string == "Lapack") {
      CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES = ON;
    } else if ( value_string == "Internal") {
      CENO_Execution_Mode::USE_LAPACK_LEAST_SQUARES = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "CENO_Additional_Approximate_Constraints") == 0) {
    value >> value_string;
    if (value_string == "Yes") {
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = ON;
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = OFF;
    } else if (value_string == "No") {
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = OFF;
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;
    }
    i_command = 0;

  } else if (strcmp(code, "Accuracy_Assessment_Near_Curved_Boundaries") == 0) {
    value >> value_string;
    if ( value_string == "Ignore_Cells") {
      CENO_Execution_Mode::IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = ON;
    } else if ( value_string == "Include_Cells") {
      CENO_Execution_Mode::IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "AMR_High_Order_Criteria") == 0) {
    value >> value_string;
    if (value_string == "Smoothness_Indicator") {
      CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = ON;
    } else {
      CENO_Execution_Mode::USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "High_Order_Message_Passing") == 0) {
    value >> value_string;
    if ( value_string == "Yes") {
      CENO_Execution_Mode::HIGH_ORDER_MESSAGE_PASSING = ON;
    } else {
      CENO_Execution_Mode::HIGH_ORDER_MESSAGE_PASSING = OFF;
    }
    i_command = 0;

  } else if (strcmp(code, "CENO_Verbose") == 0) {
    value >> value_string;
    if (value_string == "Yes" || value_string == "YES") {
      CENO_Execution_Mode::CENO_VERBOSE = ON;
      i_command = 0;
    } else if ( value_string == "No" || value_string == "NO") {
      CENO_Execution_Mode::CENO_VERBOSE = OFF;
      i_command = 0;
    } else {
      i_command = INVALID_INPUT_VALUE;
    }

  /* ----- Parse input to set the static variables in the CENO_Tolerances class ----- */

  } else if (strcmp(code, "CENO_Epsilon") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::epsilon;
  } else if (strcmp(code, "CENO_Absolute_Epsilon") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::epsilon_absolute;
  } else if (strcmp(code, "CENO_Relative_Epsilon") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::epsilon_relative;
  } else if (strcmp(code, "CENO_Tolerance") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::Fit_Tolerance;
  } else if (strcmp(code, "CENO_NonSensitivity") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::Fit_Tolerance_NonSensitivity;
  } else if (strcmp(code, "CENO_AMR_Units") == 0) {
    i_command = 0;
    value >> CENO_Tolerances::AMR_Smoothness_Units;
  } else {
    i_command = INVALID_INPUT_CODE;
    cout << "\n Default CENO_Execution_Mode and CENO_Tolerance paramaters in use. \n";
    cout.flush();
  } // endif

  // Update all dependent tolerances
  CENO_Tolerances::UpdateDependentTolerances();

  return i_command;

}
