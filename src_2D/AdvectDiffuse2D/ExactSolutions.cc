/*!\file SourceFields.cc
  \brief Source file initializing/implementing member variables/functions that belong to classes defined in SourceFields.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "AdvectDiffuse2DInput.h"
#include "AdvectDiffuse2DExactSolutions.h"


ExactSolutionBasicType::ExactSolutionBasicType(void): ExactSolutionName("Not named"), Accuracy_Parameter(Soln) {  };
ExactSolutionBasicType::~ExactSolutionBasicType(void){  };

/************************************
 * ExactSolutionBasicType Members   *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
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
      i_command = INVALID_INPUT_CODE;
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

/************************************
 * Laplace_I_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Laplace_I_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								 int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Laplace_I_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Laplace_I_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laplace_I_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Laplace_II_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Laplace_II_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								  int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Laplace_II_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Laplace_II_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laplace_II_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Laplace_III_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Laplace_III_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								   int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Laplace_III_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Laplace_III_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laplace_III_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Laplace_IV_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Laplace_IV_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								  int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "mu_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> mu;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Laplace_IV_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> mu : " << mu;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Laplace_IV_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laplace_IV_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&mu,
			1, 
			MPI::DOUBLE, 0);
 
#endif
}

/************************************
 * Laplace_V_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Laplace_V_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								 int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "D_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> D;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "mu_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> mu;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Laplace_V_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> D : " << D
	   << "\n     -> mu : " << mu;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Laplace_V_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laplace_V_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&D,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&mu,
			1, 
			MPI::DOUBLE, 0);
 
#endif
}

/************************************
 * Poisson_I_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Poisson_I_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								 int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "a_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "beta_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Poisson_I_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> a : " << a
	   << "\n     -> beta : " << beta;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Poisson_I_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Poisson_I_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Poisson_II_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Poisson_II_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								  int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "a_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "beta_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif
}

/*! 
 * Print relevant parameters
 */
void Poisson_II_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> a : " << a
	   << "\n     -> beta : " << beta;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Poisson_II_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Poisson_II_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Poisson_III_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Poisson_III_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								   int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "a_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "beta_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Poisson_III_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> a : " << a
	   << "\n     -> beta : " << beta;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Poisson_III_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Poisson_III_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/************************************
 * Poisson_IV_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Poisson_IV_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								  int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "a_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "beta_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Poisson_IV_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> a : " << a
	   << "\n     -> beta : " << beta;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Poisson_IV_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Poisson_IV_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1, 
			MPI::DOUBLE, 0);
 
#endif
}

/************************************
 * Poisson_V_ExactSolution Members  *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void Poisson_V_ExactSolution::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								 int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "A_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> A;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "B_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> B;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "C_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "a_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> a;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "beta_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> C;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else {
    i_command = INVALID_INPUT_CODE;
    return;
  } // endif

}

/*! 
 * Print relevant parameters
 */
void Poisson_V_ExactSolution::Print_Info(std::ostream & out_file){
  out_file << "\n     -> A : " << A
	   << "\n     -> B : " << B
	   << "\n     -> C : " << C
	   << "\n     -> a : " << a
	   << "\n     -> beta : " << beta;

  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the Poisson_V_ExactSolution variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Poisson_V_ExactSolution::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&A,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&B,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&C,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&a,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&beta,
			1, 
			MPI::DOUBLE, 0);
 
#endif
}


/************************************************************
 * StationaryHeatEqnWithLinearSource_ExactSolution Members  *
 ***********************************************************/

/*! 
 * Parse next input control parameter
 */
void StationaryHeatEqnWithLinearSource_ExactSolution::
Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
				   int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "lambda_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> lambda;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CoordA_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CoordA;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "SolnA_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> SolnA;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "CoordB_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CoordB;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "SolnB_Coeff") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> SolnB;
    IP.Input_File.getline(buffer, sizeof(buffer));

  } else if (strcmp(IP.Next_Control_Parameter, "Direction_Of_Variation") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "x-axis") == 0 ){
      Direction_Of_Variation = X_DIRECTION;
    } else if ( strcmp(IP.Next_Control_Parameter, "y-axis") == 0 ) {
      Direction_Of_Variation = Y_DIRECTION;
    } else {
      i_command = INVALID_INPUT_CODE;
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
void StationaryHeatEqnWithLinearSource_ExactSolution::
Print_Info(std::ostream & out_file){
  out_file << "\n     -> lambda : " << lambda;
  if (Direction_Of_Variation == X_DIRECTION){
    out_file << "\n     -> Dirichlet BC in x-direction"
	     << "\n         -> X-CoordA: " << CoordA
	     << "\n         -> Soln(X-CoordA): " << SolnA
	     << "\n         -> X-CoordB: " << CoordB
	     << "\n         -> Soln(X-CoordB): " << SolnB;
  } else {
    out_file << "\n     -> Dirichlet BC in y-direction"
	     << "\n         -> Y-CoordA: " << CoordA
	     << "\n         -> Soln(Y-CoordA): " << SolnA
	     << "\n         -> Y-CoordB: " << CoordB
	     << "\n         -> Soln(Y-CoordB): " << SolnB;
  }
  
  // call the base Print_Info
  ExactSolutionBasicType::Print_Info(out_file);
}

/*!
 * Broadcast the StationaryHeatEqnWithLinearSource_ExactSolution
 * variables to all processors associated with the specified 
 * communicator from the specified processor using the MPI broadcast 
 * routine.
 */
void StationaryHeatEqnWithLinearSource_ExactSolution::
Broadcast(void){

#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&lambda,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&CoordA,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&SolnA,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&CoordB,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&SolnB,
			1, 
			MPI::DOUBLE, 0);

  // Update all parameters on processors different than
  // the main one
  if (!CFFC_Primary_MPI_Processor()){
    Set_ParticularSolution_Parameters();
  }
 
#endif
}
