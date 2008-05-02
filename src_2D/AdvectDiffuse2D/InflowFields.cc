/*!\file InflowFields.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in InflowFields.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "AdvectDiffuse2DInput.h"
#include "InflowFields.h"


InflowFieldBasicType::InflowFieldBasicType(void): InflowFieldName("Not named") {  };
InflowFieldBasicType::~InflowFieldBasicType(void){  };


/**************************************
 * Sinusoidal_II_InflowField Members  *
 *************************************/

/*! 
 * Parse next input control parameter
 */
void Sinusoidal_II_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								   int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Reference_Point") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferencePoint;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Sinusoidal_II_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> Reference Point : " << ReferencePoint;
}

/*!
 * Broadcast the Sinusoidal_II_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Sinusoidal_II_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
			1, 
			MPI::DOUBLE, 0);

#endif
}



/**************************************
 * Sinusoidal_III_InflowField Members  *
 *************************************/

/*! 
 * Parse next input control parameter
 */
void Sinusoidal_III_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								    int & i_command){
  
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Reference_Point") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferencePoint;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_R_Min_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> R_Min;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_R_Max_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> R_Max;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Sinusoidal_III_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> Reference Point : " << ReferencePoint
	   << "\n     -> R_Min : " << R_Min
	   << "\n     -> R_Max : " << R_Max;
}

/*!
 * Broadcast the Sinusoidal_III_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Sinusoidal_III_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&R_Min,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&R_Max,
			1, 
			MPI::DOUBLE, 0);

#endif
}


/**************************************
 * Sinusoidal_IV_InflowField Members  *
 *************************************/

/*! 
 * Parse next input control parameter
 */
void Sinusoidal_IV_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								   int & i_command){
  
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Reference_Point") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferencePoint;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_R_Min_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> R_Min;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_R_Max_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> R_Max;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Sinusoidal_IV_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> Reference Point : " << ReferencePoint
	   << "\n     -> R_Min : " << R_Min
	   << "\n     -> R_Max : " << R_Max;
}

/*!
 * Broadcast the Sinusoidal_IV_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Sinusoidal_IV_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&R_Min,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&R_Max,
			1, 
			MPI::DOUBLE, 0);

#endif
}


/*********************************
 * Constant_InflowField Members  *
 ********************************/

/*! 
 * Parse next input control parameter
 */
void Constant_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
							      int & i_command){
  
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];
  
  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Constant_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A;
}

/*!
 * Broadcast the Constant_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Constant_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
#endif
}

/*********************************************
 * Hyperbolic_Tangent_I_InflowField Members  *
 ********************************************/

/*! 
 * Parse next input control parameter
 */
void Hyperbolic_Tangent_I_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
									  int & i_command){
  
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Reference_Point") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferencePoint;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_Magnitude_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Magnitude;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_Steepness_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Steepness;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Hyperbolic_Tangent_I_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> M : " << Magnitude
	   << "\n     -> S : " << Steepness
	   << "\n     -> r0: " << ReferencePoint;
}

/*!
 * Broadcast the Hyperbolic_Tangent_I_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Hyperbolic_Tangent_I_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Magnitude,
			1, 
			MPI::DOUBLE, 0);

  MPI::COMM_WORLD.Bcast(&Steepness,
			1, 
			MPI::DOUBLE, 0);
#endif
}

/*************************************************************
 * Squared_Exponential_Times_Sinusoidal_InflowField Members  *
 ************************************************************/

/*! 
 * Parse next input control parameter
 */
void Squared_Exponential_Times_Sinusoidal_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
											  int & i_command){
  
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Reference_Point") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> ReferencePoint;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_Magnitude_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Magnitude;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_Steepness_Frequency_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> SteepnessFrequency;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Inflow_Power_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> Power;
    IP.Input_File.getline(buffer, sizeof(buffer));
    if (Power < 0) i_command = INVALID_INPUT_VALUE;

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Squared_Exponential_Times_Sinusoidal_InflowField::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> M : " << Magnitude
	   << "\n     -> S : " << SteepnessFrequency
	   << "\n     -> P : " << Power
	   << "\n     -> r0: " << ReferencePoint;
}

/*!
 * Broadcast the Squared_Exponential_Times_Sinusoidal_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Squared_Exponential_Times_Sinusoidal_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Magnitude,
			1, 
			MPI::DOUBLE, 0);

  MPI::COMM_WORLD.Bcast(&SteepnessFrequency,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Power,
			1, 
			MPI::INT, 0);
#endif
}
