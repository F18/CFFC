/*!\file SourceFields.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in SourceFields.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "SourceFields.h"
#include "AdvectDiffuse2DInput.h"

SourceFieldBasicType::SourceFieldBasicType(void): FieldName("Not named"){  }
SourceFieldBasicType::~SourceFieldBasicType(void){  }
SourceFieldThatRequireIntegration::~SourceFieldThatRequireIntegration(void){ }
SourceFieldThatDoesNotRequireIntegration::~SourceFieldThatDoesNotRequireIntegration(void){ }


/*******************************
 * Linear_SourceField Members  *
 ******************************/

/*! 
 * Parse next input control parameter
 */
void Linear_SourceField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
							    int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Source_Linear_Tau_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> tau;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Linear_SourceField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> Relaxation Time : " << tau;
}

/*!
 * Broadcast the Linear_SourceField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Linear_SourceField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&tau,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}


/************************************
 * Exponential_SourceField Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Exponential_SourceField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								 int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Source_Exponential_A") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Source_Exponential_Beta") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> beta;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Exponential_SourceField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> Linear Coefficient : " << a;
  out_file << "\n     -> Exponent Coefficient : " << beta;
}

/*!
 * Broadcast the Exponential_SourceField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Exponential_SourceField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1,
			MPI::DOUBLE, 0);
 
#endif
}
