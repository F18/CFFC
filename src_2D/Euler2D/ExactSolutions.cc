/*!\file ExactSolutions.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in ExactSolutions.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "Euler2DInput.h"
#include "ExactSolutions.h"


ExactSolutionBasicType::ExactSolutionBasicType(void): ExactSolutionName("Not named"), Accuracy_Parameter(Soln) {  };
ExactSolutionBasicType::~ExactSolutionBasicType(void){  };

/************************************
 * ExactSolutionBasicType Members   *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP,
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
void ExactSolutionBasicType::Print_Info(std::ostream & out_file){
  if (Accuracy_Parameter == Soln){
    out_file << "\n     -> Measurement of Accuracy Based on Exact Solution";
  } else {
    out_file << "\n     -> Measurement of Accuracy Based on Exact Gradient";
  }
}

/*!
 * Broadcast the ExactSolutionBasicType variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void ExactSolutionBasicType::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&Accuracy_Parameter,
			1, 
			MPI::INT, 0);

#endif
}

/***************************************
 * Ringleb_Flow_ExactSolution Members  *
 **************************************/

/*! 
 * Parse next input control parameter
 */
void Ringleb_Flow_ExactSolution::Parse_Next_Input_Control_Parameter(Euler2D_Input_Parameters & IP,
								    int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
#if 0             		// Replace "..." with the proper variables
  if (strcmp(IP.Next_Control_Parameter, "...") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ...;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "...") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ...;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
#endif
}

/*! 
 * Print relevant parameters
 */
void Ringleb_Flow_ExactSolution::Print_Info(std::ostream & out_file){
#if 0             		// Replace "..." with the proper variables
  out_file << "\n     -> "..." : " << "..."
	   << "\n     -> "..." : " << "...";

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
#endif
}

/*!
 * Broadcast the Ringleb_Flow_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Ringleb_Flow_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
#if 0             		// Replace "..." with the proper variables
  MPI::COMM_WORLD.Bcast(&"...",
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&"...",
			1, 
			MPI::DOUBLE, 0);
#endif
  
#endif
}

