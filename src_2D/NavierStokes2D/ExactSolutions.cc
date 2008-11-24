/*!\file ExactSolutions.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in ExactSolutions.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "NavierStokes2DInput.h"
#include "ExactSolutions.h"

ExactSolutionBasicType::ExactSolutionBasicType(void): ExactSolutionName("Not named"), Accuracy_Parameter(Soln) {  };
ExactSolutionBasicType::~ExactSolutionBasicType(void){  };

/************************************
 * ExactSolutionBasicType Members   *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
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


/*******************************************
 * Abgrall_Function_ExactSolution Members  *
 ******************************************/

/*! 
 * Parse next input control parameter
 */
void Abgrall_Function_ExactSolution::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }  
}

/*! 
 * Print relevant parameters
 */
void Abgrall_Function_ExactSolution::Print_Info(std::ostream & out_file){

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Abgrall_Function_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Abgrall_Function_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  // None
#endif
}

/**********************************************
 * Sinusoidal_Function_ExactSolution Members  *
 *********************************************/

/*! 
 * Parse next input control parameter
 */
void Sinusoidal_Function_ExactSolution::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									   int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Min_Domain_Limit") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> DomainMinLimit;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Max_Domain_Limit") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> DomainMaxLimit;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Reference_Value") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferenceValue;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Velocity;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Direction") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Xdir") == 0 ){
      Direction = X_DIRECTION;
      i_command = 0;      
    } else if ( strcmp(IP.Next_Control_Parameter, "Ydir") == 0 ){
      Direction = Y_DIRECTION;
      i_command = 0;
    } else {
      i_command = INVALID_INPUT_CODE;
    }

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Sinusoidal_Function_ExactSolution::Print_Info(std::ostream & out_file){

  if (Direction == X_DIRECTION){
    out_file << "\n     -> Direction : X";
  } else {
    out_file << "\n     -> Direction : Y";
  }

  out_file << "\n     -> Min_Domain_Limit : " << DomainMinLimit
	   << "\n     -> Max_Domain_Limit : " << DomainMaxLimit
	   << "\n     -> Reference Value : " << ReferenceValue
	   << "\n     -> Velocity : " << Velocity;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);

}

/*!
 * Broadcast the Sinusoidal_Function_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Sinusoidal_Function_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&DomainMinLimit,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&DomainMaxLimit,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferenceValue,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Velocity,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Direction,
			1, 
			MPI::SHORT, 0);
  
#endif
}


/******************************************
 * CosSin_Function_ExactSolution Members  *
 *****************************************/

/*! 
 * Parse next input control parameter
 */
void CosSin_Function_ExactSolution::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
								       int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Min_Domain_Limit") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> DomainMinLimit;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Max_Domain_Limit") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> DomainMaxLimit;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Reference_Value") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferenceValue;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Amplitude") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Amplitude;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Velocity;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Direction") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Xdir") == 0 ){
      Direction = X_DIRECTION;
      i_command = 0;      
    } else if ( strcmp(IP.Next_Control_Parameter, "Ydir") == 0 ){
      Direction = Y_DIRECTION;
      i_command = 0;
    } else {
      i_command = INVALID_INPUT_CODE;
    }

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void CosSin_Function_ExactSolution::Print_Info(std::ostream & out_file){

  if (Direction == X_DIRECTION){
    out_file << "\n     -> Direction : X";
  } else {
    out_file << "\n     -> Direction : Y";
  }

  out_file << "\n     -> Min_Domain_Limit : " << DomainMinLimit
	   << "\n     -> Max_Domain_Limit : " << DomainMaxLimit
	   << "\n     -> Reference Value : " << ReferenceValue
	   << "\n     -> Amplitude : " << Amplitude
	   << "\n     -> Velocity : " << Velocity;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);

}

/*!
 * Broadcast the CosSin_Function_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void CosSin_Function_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&DomainMinLimit,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&DomainMaxLimit,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferenceValue,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Amplitude,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Velocity,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Direction,
			1, 
			MPI::SHORT, 0);
  
#endif
}


/************************************************
 * UnitTest_Function_ExactSolution Members  *
 ***********************************************/

/*! 
 * Parse next input control parameter
 */
void UnitTest_Function_ExactSolution::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									 int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

}

/*! 
 * Print relevant parameters
 */
void UnitTest_Function_ExactSolution::Print_Info(std::ostream & out_file){

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);

}

/*!
 * Broadcast the UnitTest_Function_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void UnitTest_Function_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
#endif
}
