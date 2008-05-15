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

  /*! This flag controls whether the CENO algorithm is used in the computation.\n
      Turn ON if you want to use the high-order CENO scheme. \n
      Turn OFF if you don't want to use this scheme. (default) \n
      ATTENTION: Some settings/memory allocations are done only if this flag is ON. \n
      ----------------------------------------------------------------------------------------  */
  static short USE_CENO_ALGORITHM;

  
  /*! Store the pseudo-inverse of the LHS term in the CENO reconstruction for every computational cell.\n
      Turn ON if you want to run in speed efficient mode. However, the memory requirements will
      increase considerably!!! (default) \n
      Turn OFF for running in the memory efficient mode. \n
      ---------------------------------------------------------------------------------------- */
  static short CENO_SPEED_EFFICIENT;


  /*! This flag controls whether the order is dropped or not in regions detected as non-smooth.\n
      Turn ON if you want the CENO scheme to produce monotone solutions in non-smooth regions. (default) \n
      Turn OFF if you want to use high-order interpolant throughout the whole computational domain. \n
      ATTENTION: If the flag is OFF the solution might blow up !!! \n
      ----------------------------------------------------------------------------------------  */
  static short CENO_DROP_ORDER;  


  /*! Use CENO smoothness indicator with padding. \n
      Turn ON to have the order dropped in all the cells surrounding a cell with non-smooth reconstruction. \n
      Turn OFF for the opposite behaviour (default) \n
      ---------------------------------------------------------------------------------------- */
  static short CENO_PADDING;

  /*! Turn ON this flag if the geometric weighting is applied. (default) \n
      Turn OFF if the geometric weighting is not applied. \n
      ---------------------------------------------------------------------------------------- */
  static short CENO_APPLY_GEOMETRIC_WEIGHTING;

  /*! Turn ON this flag if the geometric weighting used in the least-squares problem is 1.0/(Distance^2). \n
      Turn OFF if the geometric weighting is 1.0/fabs(Distance) (default) \n
      ---------------------------------------------------------------------------------------- */
  static short CENO_SQUARE_GEOM_WEIGHTING;


  /*! Turn ON if you want to include the geometric weights used in the reconstruction process
      as part of the computation of the smoothness indicator. \n
      This works only if the CENO_SUPER_FAST is turned ON (because the geometric weights are saved). \n
      Turn OFF for treating all the terms in the residual sum (SS_Residual) and regression sum (SS_Regression)
      equally (i.e. each of them has a geoemtric weight equals to 1 in the computation of the 
      smoothness indicator). (default) \n
     ---------------------------------------------------------------------------------------- */
  static short CENO_CONSIDER_WEIGHTS;

  /*! Use the piecewise constant state(s) to solve the Riemann problem at the interface of two cells 
      if any/both of the interface states calculated with the reconstruction interpolant
      has/have non-physical values (e.g. negative pressure or density). \n
      Turn ON if this feature is desired. (default) \n
      Turn OFF if you don't want to use this feature. \n
      --------------------------------------------------------------------------------------- */
  static short FORCE_WITH_PIECEWISE_CONSTANT_AT_INTERFACE;

  /*! This flag is used to switch between two methods of performing the CENO reconstruction. \n
      The first method tries to reduce the number of necessary ghost cells by using message
      passing to the overlapping layers of ghost cells in order to decide the smoothness properties
      of the cells used to compute the interface flux. This method has not been explored in practice yet. \n
      The second method uses enough layers of ghost cells to compute all necessary reconstructions
      for either flux calculation or smoothness properties assessment. \n
      Turn ON if the message passing is going to be used (i.e. first method). Not implemented yet!!! \n
      Turn OFF if no message passing is used in the reconstruction process. (default) \n
      \note Message passing might still be used for cell reconstructions at boundaries between
      blocks with mesh resolution change! This is subject to a different flag. \n
      --------------------------------------------------------------------------------------- */
  static short CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING;

  /*! Use the only the first neighbour cells to calculate the value of the smoothness indicator,
      instead of using all the cells used in the reconstruction process. \n
      Turn ON if this feature is desired. \n
      Turn OFF if you don't want to use this feature.  (default) \n
      --------------------------------------------------------------------------------------- */
  static short CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS;

  /*! Generate approximate constraints from the neighbour cells that are near constrained boundaries. \n
      Turn ON if this feature is desired. (default) \n
      Turn OFF if you don't want to use this feature. \n
      --------------------------------------------------------------------------------------- */
  static short CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS;

  /*! Extend the stencil to get more equations instead of using approximate constraints with
      the boundary conditions of the neighbour cells that are near constrained boundaries. \n
      Turn ON if this feature is desired. \n
      Turn OFF if you don't want to use this feature. (default) \n
      --------------------------------------------------------------------------------------- */
  static short CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL;

  /*! Error functions might be hard or impossible to be integrated accurately over domains with
      curved boundaries. Therefore, this flag controls whether the cells near these boundaries 
      are included or ignored. \n
      Turn ON if this feature is desired. \n
      Turn OFF if you don't want to use this feature. (default) \n
      --------------------------------------------------------------------------------------- */  
  static short IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT;

  /*! Error functions might be hard or impossible to be integrated accurately over domains with
      curved boundaries. Therefore, this flag controls whether the cells near these boundaries 
      are included or ignored. \n
      Turn ON if this feature is desired. \n
      Turn OFF if you don't want to use this feature. (default) \n
      --------------------------------------------------------------------------------------- */  
  static short USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA;  

  /*! The high-order k-Exact reconstruction involves the solution of a linear least-squares problem,
      which can be obtained with different subroutines. Currently, a Lapack subroutine and an 
      internal one (L. Ivan's implementation) can be used. The Lapack one is faster but assumes the 
      rank of the matrix equal to the number of columns whereas the internal one can be configured to
      detect ranks smaller than this (i.e. it still gives a solution) but it's slower. \n
      Turn ON to use the Lapack library subroutine. (default) \n
      Turn OFF to use the internal one.  
      --------------------------------------------------------------------------------------- */
  static short USE_LAPACK_LEAST_SQUARES;

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
      i_command = INVALID_INPUT_VALUE;
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
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Apply_Geom_Weighting") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_APPLY_GEOMETRIC_WEIGHTING = ON;
    } else {
      CENO_APPLY_GEOMETRIC_WEIGHTING = OFF;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Use_Geom_Weighting_For_Smoothness_Analysis") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_CONSIDER_WEIGHTS = ON;
    } else {
      CENO_CONSIDER_WEIGHTS = OFF;
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
  } else if (strcmp(IP.Next_Control_Parameter, "CENO_With_Message_Passing") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = ON;
      throw runtime_error("CENO_Execution_Mode ERROR! CENO reconstruction with message passing has not been implemented yet!");
    } else {
      CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING = OFF;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Smoothness_Indicator") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Use_First_Neighbours") == 0 ){
      CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = ON;
    } else if ( strcmp(IP.Next_Control_Parameter, "Use_All_Neighbours") == 0 ){
      CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS = OFF;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Least_Squares_Solver") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Lapack") == 0 ){
      USE_LAPACK_LEAST_SQUARES = ON;
    } else if ( strcmp(IP.Next_Control_Parameter, "Internal") == 0 ){
      USE_LAPACK_LEAST_SQUARES = OFF;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "CENO_Additional_Approximate_Constraints") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Yes") == 0 ){
      CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = ON;
      CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = OFF;
    } else if ( strcmp(IP.Next_Control_Parameter, "No") == 0 ){
      CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS = OFF;
      CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL = ON;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "Accuracy_Assessment_Near_Curved_Boundaries") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Ignore_Cells") == 0 ){
      IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = ON;
    } else if ( strcmp(IP.Next_Control_Parameter, "Include_Cells") == 0 ){
      IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT = OFF;
    }
    i_command = 0;

  } else if (strcmp(IP.Next_Control_Parameter, "AMR_High_Order_Criteria") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Smoothness_Indicator") == 0 ){
      USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = ON;
    } else {
      USE_SMOOTHNESS_INDICATOR_FOR_AMR_CRITERIA = OFF;
    }
    i_command = 0;

  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}

#endif
