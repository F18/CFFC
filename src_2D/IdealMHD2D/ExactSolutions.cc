/*!\file ExactSolutions.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in ExactSolutions.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "MHD2DInput.h"
#include "ExactSolutions.h"

ExactSolutionBasicType_MHD2D::ExactSolutionBasicType_MHD2D(void): ExactSolutionName("Not named"),
								  Accuracy_Parameter(Soln) {  };
ExactSolutionBasicType_MHD2D::~ExactSolutionBasicType_MHD2D(void){  };

/********************************************
 * ExactSolutionBasicType_MHD2D Members   *
 *******************************************/

/*! 
 * Parse next input control parameter
 */
void ExactSolutionBasicType_MHD2D::Parse_Next_Input_Control_Parameter(MHD2D_Input_Parameters & IP,
								      int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];
  
  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Base_Accuracy_Measurement_On") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Solution") == 0 ){
      Accuracy_Parameter = Soln;
    } else if ( strcmp(IP.Next_Control_Parameter, "Gradient") == 0 ) {
      Accuracy_Parameter = Grad;
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void ExactSolutionBasicType_MHD2D::Print_Info(std::ostream & out_file){
  if (Accuracy_Parameter == Soln){
    out_file << "\n     -> Measurement of Accuracy Based on Exact Solution";
  } else {
    out_file << "\n     -> Measurement of Accuracy Based on Exact Gradient";
  }
}

/*!
 * Broadcast the ExactSolutionBasicType_MHD2D variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void ExactSolutionBasicType_MHD2D::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&Accuracy_Parameter,
			1, 
			MPI::INT, 0);

#endif
}


/************************************************
 * UnitTest_Function_ExactSolution_MHD Members  *
 ***********************************************/

/*! 
 * Parse next input control parameter
 */
void UnitTest_Function_ExactSolution_MHD::Parse_Next_Input_Control_Parameter(MHD2D_Input_Parameters & IP,
									     int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_MHD2D::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];
  
}

/*! 
 * Print relevant parameters
 */
void UnitTest_Function_ExactSolution_MHD::Print_Info(std::ostream & out_file){
  
  // call the base Print_Info
  ExactSolutionBasicType_MHD2D::Print_Info(out_file);
  
}

/*!
 * Broadcast the UnitTest_Function_ExactSolution_MHD variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void UnitTest_Function_ExactSolution_MHD::Broadcast(void){
#ifdef _MPI_VERSION
  
#endif
}
