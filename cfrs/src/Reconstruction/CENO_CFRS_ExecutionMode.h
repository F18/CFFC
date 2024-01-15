/*!\file CENO_ExecutionMode.h
  \brief Definition of flags that control the execution of CENO high-order reconstruction. */

#ifndef _CENO_CFRS_EXECUTIONMODE_INCLUDED
#define _CENO_CFRS_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../../../src_2D/Utilities/Utilities.h"
#include "../../../src_2D/CFD/CFD.h"

/*!
 * Class: CENO_CFRS_Execution_Mode
 *
 * @brief Definition of flags that control the execution of CENO high-order reconstruction
 *
 */
class CENO_CFRS_Execution_Mode{
  
public:

  // set all flags to default values
  static void SetDefaults(void);

   
  /*! Store the pseudo-inverse of the LHS term in the CENO reconstruction for every computational cell.\n
      Turn ON if you want to run in speed efficient mode. However, the memory requirements will
      increase considerably!!! (default) \n
      Turn OFF for running in the memory efficient mode. \n
      ---------------------------------------------------------------------------------------- */
  static short USE_PSEUDO_INVERSE;

 /*! Reduce the order of reconstruction to limited piecewise linear reconstruction in the cells 
     flagged as under-resolved or non-smooth by the solution smoothness indicator. \n
     Turn ON if you want to recompute the reconstruction in the flagged cells using limited 
     piecewise linear reconstruction. (default) \n
     Turn OFF if you want to keep the unlimited high-order reconstruction. \n
   ---------------------------------------------------------------------------------------- */
  static short REDUCE_ORDER;


  /*! Turn ON this flag if the geometric weighting is applied. (default) \n
      Turn OFF if the geometric weighting is not applied. \n
      ---------------------------------------------------------------------------------------- */
  //  static short CENO_APPLY_GEOMETRIC_WEIGHTING;
  
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  // template<class Input_Parameters_Type>
  // static void Get_Next_Input_Control_Parameter(Input_Parameters_Type &IP);

  static void Print_Info(std::ostream & out_file);

  //  static void Broadcast(void);

protected:
  CENO_CFRS_Execution_Mode(void);   //!< Private default constructor
  CENO_CFRS_Execution_Mode(const CENO_CFRS_Execution_Mode&); //!< Private copy constructor
  CENO_CFRS_Execution_Mode& operator=(const CENO_CFRS_Execution_Mode&); //!< Private assignment operator

};

//! Parse the input control parameters for 
//  settings related to CENO_CFRS_Execution_Mode class
template<class Input_Parameters_Type> inline
void CENO_CFRS_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
    if (i_command != INVALID_INPUT_CODE){
      return;
    }

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Pseudo_Inverse") == 0){
    i_command = 36;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Pseudo_Inverse, 
	   IP.Next_Control_Parameter);
    
    // Add the parse word for the new function //
    if (strcmp(IP.Pseudo_Inverse, "Yes") == 0)
      CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE = ON;
    else if (strcmp(IP.Pseudo_Inverse, "No") == 0)
      CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE = OFF;
    else if (strcmp(IP.Pseudo_Inverse, "ON") == 0)
      CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE = ON;
    else if (strcmp(IP.Pseudo_Inverse, "OFF") == 0)
      CENO_CFRS_Execution_Mode::USE_PSEUDO_INVERSE = OFF;

  } else if (strcmp(IP.Next_Control_Parameter, "Reduce_Order") == 0){
    i_command = 37;
    Get_Next_Input_Control_Parameter(IP);
    strcpy(IP.Reduce_Order, 
	   IP.Next_Control_Parameter);
    
    // Add the parse word for the new function //
    if (strcmp(IP.Reduce_Order, "Yes") == 0)
      CENO_CFRS_Execution_Mode::REDUCE_ORDER = ON;
    else if (strcmp(IP.Reduce_Order, "No") == 0)
      CENO_CFRS_Execution_Mode::REDUCE_ORDER = OFF;
    else if (strcmp(IP.Reduce_Order, "ON") == 0)
      CENO_CFRS_Execution_Mode::REDUCE_ORDER = ON;
    else if (strcmp(IP.Reduce_Order, "OFF") == 0)
      CENO_CFRS_Execution_Mode::REDUCE_ORDER = OFF;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

/********************************************************
 * Routine: Get_Next_Input_Control_Parameter            *
 *                                                      *
 * Get the next input control parameter from the input  *
 * file.                                                *
 *                                                      *
 ********************************************************/
//template<class Input_Parameters_Type> inline
//void CENO_CFRS_Execution_Mode::Get_Next_Input_Control_Parameter(Input_Parameters_Type &IP) {
//
//    int i;
//    char buffer[256];
//
//    IP.Line_Number = IP.Line_Number + 1;
//    IP.Input_File.getline(buffer, sizeof(buffer));
//    i = 0;
//    if (buffer[0] != '#') {
//       while (1) {
//          if (buffer[i] == ' ' || buffer[i] == '=' ) break;
//          i = i + 1;
//          if (i > strlen(buffer) ) break;
//       } /* endwhile */
//       buffer[i] = '\0';
//    } /* endif */
//    strcpy(IP.Next_Control_Parameter, buffer);
//}


#endif
