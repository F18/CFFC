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

ExactSolutionBasicType_NavierStokes2D::ExactSolutionBasicType_NavierStokes2D(void): ExactSolutionName("Not named"),
										    Accuracy_Parameter(Soln) {  };
ExactSolutionBasicType_NavierStokes2D::~ExactSolutionBasicType_NavierStokes2D(void){  };

/************************************
 * ExactSolutionBasicType_NavierStokes2D Members   *
 ***********************************/

/*! 
 * Parse next input control parameter
 */
void ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
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
void ExactSolutionBasicType_NavierStokes2D::Print_Info(std::ostream & out_file){
  if (Accuracy_Parameter == Soln){
    out_file << "\n     -> Measurement of Accuracy Based on Exact Solution";
  } else {
    out_file << "\n     -> Measurement of Accuracy Based on Exact Gradient";
  }
}

/*!
 * Broadcast the ExactSolutionBasicType_NavierStokes2D variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void ExactSolutionBasicType_NavierStokes2D::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&Accuracy_Parameter,
			1, 
			MPI::INT, 0);

#endif
}

/*! Output the variable names that defined the exact solution in a format 
 *  suitable for Tecplot to the provided output stream.
 *  The default output assumes that only four variables (i.e. rho, u, v, p)
 *  are known exactly.
 */
void ExactSolutionBasicType_NavierStokes2D::Output_Tecplot_Title(std::ostream & out_file) const{
  out_file << "\"ExactSoln_rho\" \\ \n"
	   << "\"ExactSoln_u\" \\ \n"
	   << "\"ExactSoln_v\" \\ \n"
	   << "\"ExactSoln_p\" \\ \n";
}

/*! Write the Tecplot double precision format for the default output. */
void ExactSolutionBasicType_NavierStokes2D::Output_Tecplot_Double_Precision(std::ostream & out_file) const{
  out_file << "DOUBLE DOUBLE DOUBLE DOUBLE ";
}

/*******************************************
 * Abgrall_Function_ExactSolution_NS Members  *
 ******************************************/

/*! 
 * Parse next input control parameter
 */
void Abgrall_Function_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }  
}

/*! 
 * Print relevant parameters
 */
void Abgrall_Function_ExactSolution_NS::Print_Info(std::ostream & out_file){

  // call the base Print_Info
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);
}

/*!
 * Broadcast the Abgrall_Function_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Abgrall_Function_ExactSolution_NS::Broadcast(void){
#ifdef _MPI_VERSION
  // None
#endif
}

/**********************************************
 * Sinusoidal_Function_ExactSolution_NS Members  *
 *********************************************/

/*! 
 * Parse next input control parameter
 */
void Sinusoidal_Function_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									      int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

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
void Sinusoidal_Function_ExactSolution_NS::Print_Info(std::ostream & out_file){

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
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);

}

/*!
 * Broadcast the Sinusoidal_Function_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Sinusoidal_Function_ExactSolution_NS::Broadcast(void){
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
 * CosSin_Function_ExactSolution_NS Members  *
 *****************************************/

/*! 
 * Parse next input control parameter
 */
void CosSin_Function_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									  int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

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
void CosSin_Function_ExactSolution_NS::Print_Info(std::ostream & out_file){

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
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);

}

/*!
 * Broadcast the CosSin_Function_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void CosSin_Function_ExactSolution_NS::Broadcast(void){
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
 * UnitTest_Function_ExactSolution_NS Members  *
 ***********************************************/

/*! 
 * Parse next input control parameter
 */
void UnitTest_Function_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									    int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

}

/*! 
 * Print relevant parameters
 */
void UnitTest_Function_ExactSolution_NS::Print_Info(std::ostream & out_file){

  // call the base Print_Info
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);

}

/*!
 * Broadcast the UnitTest_Function_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void UnitTest_Function_ExactSolution_NS::Broadcast(void){
#ifdef _MPI_VERSION
  
#endif
}


/************************************************
 * Viscous_Channel_Flow_ExactSolution_NS Members  *
 ***********************************************/

/*! 
 * Parse next input control parameter
 */
void Viscous_Channel_Flow_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									       int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

}

/*! 
 * Print relevant parameters
 */
void Viscous_Channel_Flow_ExactSolution_NS::Print_Info(std::ostream & out_file){

  out_file << "\n     -> Channel Height : " << HeightDomain
	   << "\n     -> Channel Length : " << LengthDomain
    	   << "\n     -> Pressure Change : " << DeltaPressure
	   << "\n     -> Reference Density : " << Wo.rho
	   << "\n     -> Reference Pressure : " << Wo.p
	   << "\n     -> Upper Wall Speed : " << Vwall.x;

  // call the base Print_Info
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);

}

/*!
 * Broadcast the Viscous_Channel_Flow_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Viscous_Channel_Flow_ExactSolution_NS::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&Wo.rho,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Wo.v.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Wo.v.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Wo.p,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&DeltaPressure,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&HeightDomain,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&LengthDomain,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Vwall.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Vwall.y,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/*!
 * Set particular solution parameters based on
 * the input parameters.
 */
void Viscous_Channel_Flow_ExactSolution_NS::
Set_ParticularSolution_Parameters(const NavierStokes2D_Input_Parameters & IP){

  ostringstream ErrorMsg;

  Wo = IP.Wo;
  DeltaPressure = IP.dp;
  switch(IP.i_Grid) {
  case GRID_CARTESIAN_UNIFORM :
    HeightDomain = IP.Box_Height;
    LengthDomain = IP.Box_Width;
    break;
  case GRID_SQUARE :
    HeightDomain = IP.Box_Width;
    LengthDomain = IP.Box_Width;
    break;
  case GRID_RECTANGULAR_BOX :
    HeightDomain = IP.Box_Height;
    LengthDomain = IP.Box_Width;
    break;
  default:
    ErrorMsg << "Viscous_Channel_Flow_ExactSolution_NS::Set_ParticularSolution_Parameters() ERROR!\n"
	     << "The current mesh is not suitable for this exact solution!\n";
    throw runtime_error(ErrorMsg.str());
  }
  Vwall = IP.Vwall;
  
}


/************************************************
 * Laminar_Flat_Plate_ExactSolution_NS Members  *
 ***********************************************/

/*! 
 * Parse next input control parameter
 */
void Laminar_Flat_Plate_ExactSolution_NS::Parse_Next_Input_Control_Parameter(NavierStokes2D_Input_Parameters & IP,
									     int & i_command){

  // Call the parser from the base class
  ExactSolutionBasicType_NavierStokes2D::Parse_Next_Input_Control_Parameter(IP,i_command);

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }
  
  char buffer[256];

}

/*! 
 * Print relevant parameters
 */
void Laminar_Flat_Plate_ExactSolution_NS::Print_Info(std::ostream & out_file){

  out_file << "\n     -> Plate length : " << plate_length
	   << "\n     -> Reference Density : " << Winf.rho
	   << "\n     -> Reference Pressure : " << Winf.p;

  // call the base Print_Info
  ExactSolutionBasicType_NavierStokes2D::Print_Info(out_file);

}

/*!
 * Broadcast the Viscous_Channel_Flow_ExactSolution_NS variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void Laminar_Flat_Plate_ExactSolution_NS::Broadcast(void){
#ifdef _MPI_VERSION

  MPI::COMM_WORLD.Bcast(&Winf.rho,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Winf.v.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Winf.v.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&Winf.p,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&plate_length,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/*!
 * Set particular solution parameters based on
 * the input parameters.
 */
void Laminar_Flat_Plate_ExactSolution_NS::
Set_ParticularSolution_Parameters(const NavierStokes2D_Input_Parameters & IP){

  Winf = IP.Wo;
  plate_length = IP.Plate_Length;
}
