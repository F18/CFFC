/*!\file CENO_ExecutionMode.h
  \brief Definition of flags that control the execution of CENO high-order reconstruction. */

#ifndef _CENO_EXECUTIONMODE_INCLUDED
#define _CENO_EXECUTIONMODE_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"

/*!
 * Class: CENO_Execution_Mode
 *
 * @brief Definition of flags that control the execution of CENO high-order reconstruction
 *
 */
class CENO_Execution_Mode{
  
public:

  // set all flags to default values
  static void SetDefaults(void);

  /* This flag controls whether the CENO algorithm is used in the computation.
     Turn ON if you want to use the high-order CENO scheme.
     Turn OFF if you don't want to use this scheme. (default)
     ATTENTION: Some settings/memory allocations are done only if this flag is ON.
     ----------------------------------------------------------------------------------------  */
  static short USE_CENO_ALGORITHM;

  
  /* Store the pseudo-inverse of the LHS term in the CENO reconstruction for every computational cell.
     Turn ON if you want to run in speed efficient mode. However, the memory requirements will
     increase considerably!!! (default)
     Turn OFF for running in the memory efficient mode.
     ---------------------------------------------------------------------------------------- */
  static short CENO_SPEED_EFFICIENT;


  /* This flag controls whether the order is dropped or not in regions detected as non-smooth.
     Turn ON if you want the CENO scheme to produce monotone solutions in non-smooth regions. (default)
     Turn OFF if you want to use high-order interpolant throughout the whole computational domain.
     ATTENTION: If the flag is OFF the solution might blow up !!! 
     ----------------------------------------------------------------------------------------  */
  static short CENO_DROP_ORDER;  


  /* Use CENO smoothness indicator with padding 
     Turn ON to have the order dropped in all the cells surrounding a cell with non-smooth reconstruction
     Turn OFF for the opposite behaviour (default)
     ---------------------------------------------------------------------------------------- */
  static short CENO_PADDING;

  /* Turn ON this flag if the geometric weighting used in the least-squares problem is 1.0/(Distance^2)
     Turn OFF if the geometric weighting is 1.0/fabs(Distance) (default)
     ---------------------------------------------------------------------------------------- */
  static short CENO_SQUARE_GEOM_WEIGHTING;


  /* Turn ON if you want to include the geometric weights used in the reconstruction process
     as part of the computation of the smoothness indicator.
     This works only if the CENO_SUPER_FAST is turned ON (because the geometric weights are saved).
     Turn OFF for treating all the terms in the residual sum (SS_Residual) and regression sum (SS_Regression)
     equally (i.e. each of them has a geoemtric weight equals to 1 in the computation of the 
     smoothness indicator). (default)
     ---------------------------------------------------------------------------------------- */
  static short CENO_CONSIDER_WEIGHTS;

  /* Use the piecewise constant state(s) to solve the Riemann problem at the interface of two cells 
     if any/both of the interface states calculated with the reconstruction interpolant
     has/have non-physical values (e.g. negative pressure or density).
     Turn ON if this feature is desired. (default)
     Turn OFF if you don't want to use this feature.
     --------------------------------------------------------------------------------------- */
  static short FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE;


  static int Limiter;   //!< the limiter used for the limited linear reconstruction performed for non-smooth solutions

  
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file);

  static void Broadcast(void);

protected:
  CENO_Execution_Mode(void);   //!< Private default constructor
  CENO_Execution_Mode(const CENO_Execution_Mode&); //!< Private copy constructor
  CENO_Execution_Mode& operator=(const CENO_Execution_Mode&); //!< Private assignment operator

};

//! Parse the input control parameters for 
//  settings related to CENO_Execution_Mode class
template<class Input_Parameters_Type> inline
void CENO_Execution_Mode::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "CENO_Execution_Mode") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Speed_Efficient") == 0 ){
      CENO_SPEED_EFFICIENT = ON;
    } else if ( strcmp(IP.Next_Control_Parameter, "Memory_Efficient") == 0 ) {
      CENO_SPEED_EFFICIENT = OFF;
    } else {
      i_command = INVALID_INPUT_CODE;
      return;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Drop_Order") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_DROP_ORDER = ON;      
    } else {
      CENO_DROP_ORDER = OFF;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Padding") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_PADDING = ON;
    } else {
      CENO_PADDING = OFF;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Square_Geom_Weighting") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_SQUARE_GEOM_WEIGHTING = ON;
    } else {
      CENO_SQUARE_GEOM_WEIGHTING = OFF;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Force_With_PWC") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = ON;
    } else {
      FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE = OFF;
    }
    i_command = 0;
  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif
