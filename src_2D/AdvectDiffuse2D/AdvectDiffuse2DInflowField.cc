/*!\file AdvectDiffuse2DInflowField.cc
  \brief Source file initializing/implementing member variables/functions of class AdvectDiffuse2D_InflowField. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "AdvectDiffuse2DInflowField.h"
#include "AdvectDiffuse2DInput.h"
#include "../CFD/CFD.h"


// ==== Member variables ====
InflowFieldBasicType * AdvectDiffuse2D_InflowField::Inflow = NULL;   //!< no associated inflow field
short AdvectDiffuse2D_InflowField::i_Inflow_Field_Type = NO_INFLOW_FIELD;  //!< value for no associated inflow field

// ==== Member functions ====
/*! Default constructor is used to create the unique object of this class */
AdvectDiffuse2D_InflowField::AdvectDiffuse2D_InflowField(void){ }

/*! Destroy the inflow field object */
void AdvectDiffuse2D_InflowField::DestroyInflowFieldObject(void){
  delete Inflow;
  Inflow = NULL;		// set to no associated inflow field
  i_Inflow_Field_Type = NO_INFLOW_FIELD; // set to value for no associated inflow field
}

/*!
 * Create and return the address of the inflow field.
 * This is a singleton class (only one instance is created).\n
 * DON'T modify this function unless you are familiar with the Meyers singleton pattern.
 */
AdvectDiffuse2D_InflowField& AdvectDiffuse2D_InflowField::getInstance(void){
  static AdvectDiffuse2D_InflowField inst;
  std::atexit(DestroyInflowFieldObject);		// schedule for deallocation at exit
  return inst;
}

/*!
 * Set the inflow field based on the required field type
 */
void AdvectDiffuse2D_InflowField::SetInflowField(const short &InflowIndex){
  // deallocate the inflow field object
  DestroyInflowFieldObject();

  // create the proper exact solution and set the index accordingly
  switch (InflowIndex){
  case NO_INFLOW_FIELD:
    // Don't do anything. The values were set by DestroyInflowFieldObject() routine.
    break;
  case SINUSOIDAL_I:
    Inflow = new Sinusoidal_I_InflowField;
    break;
    
  case SINUSOIDAL_II:
    Inflow = new Sinusoidal_II_InflowField;
    break;
    
  case SINUSOIDAL_III:
    Inflow = new Sinusoidal_III_InflowField;
    break;

  case SINUSOIDAL_IV:
    Inflow = new Sinusoidal_IV_InflowField;
    break;

  case CONSTANT_INFLOW_FIELD:
    Inflow = new Constant_InflowField;
    break;

  case HYPERBOLIC_TANGENT_I:
    Inflow = new Hyperbolic_Tangent_I_InflowField;
    break;

  case EXPONENTIAL_SINUSOIDAL:
    Inflow = new Squared_Exponential_Times_Sinusoidal_InflowField;
    break;

  case TOP_HAT:
    Inflow = new Top_Hat_InflowField;
    break;
      
  case DISCONTINUOUS_EXPONENTIAL_SINUSOIDAL:
    Inflow = new Discontinuous_Exponential_Times_Sinusoidal_InflowField;
    break;

  default:
    throw runtime_error("AdvectDiffuse2D_InflowField::SetInflowField() ERROR! Unknown inflow field index.");
  }
  
  // Store exact solution type
  i_Inflow_Field_Type = InflowIndex;
}

/*!
 * Parse the input control parameters for settings 
 * related to AdvectDiffuse2D_InflowField class.
 * The first entry related to the exact solution in the
 * input file MUST specify the type of the exact solution.
 * Consequent parameters will be parsed with the parser 
 * particular to the specified inflow field type.
 */
void AdvectDiffuse2D_InflowField::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
								     int & i_command){
  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Inflow_Field") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Inflow_Sinusoidal_I") == 0 ){
      SetInflowField(SINUSOIDAL_I);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_Sinusoidal_II") == 0 ) {
      SetInflowField(SINUSOIDAL_II);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_Sinusoidal_III") == 0 ) {
      SetInflowField(SINUSOIDAL_III);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_Sinusoidal_IV") == 0 ) {
      SetInflowField(SINUSOIDAL_IV);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_Constant") == 0 ) {
      SetInflowField(CONSTANT_INFLOW_FIELD);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_HyperTan_I") == 0 ) {
      SetInflowField(HYPERBOLIC_TANGENT_I);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_ExpSinusoidal") == 0 ) {
      SetInflowField(EXPONENTIAL_SINUSOIDAL);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_TopHat") == 0 ) {
      SetInflowField(TOP_HAT);
    } else if ( strcmp(IP.Next_Control_Parameter, "Inflow_DiscontinuousExpSinusoidal") == 0 ) {
      SetInflowField(DISCONTINUOUS_EXPONENTIAL_SINUSOIDAL);
    } else {
      i_command = INVALID_INPUT_CODE;
      return;
    }
    i_command = 0;

  } else {
    // Continue parsing with the parser of the current inflow field
    if (Inflow != NULL){
      Inflow->Parse_Next_Input_Control_Parameter(IP,i_command);
    }
    return;
  } // endif

}

/*!
 * Print the relevant parameters of the AdvectDiffuse2D_InflowField class for the 
 * selected inflow field type to the provided output stream.
 */
void AdvectDiffuse2D_InflowField::Print_Info(std::ostream & out_file){

  out_file << "\n  -> Inflow field: " ; 
  if (Inflow != NULL){
    // Output field name
    out_file << Inflow->whatInflowField();

    // Output exact solution characteristic parameters
    Inflow->Print_Info(out_file);
  } else {
    // There is no associated exact solution
    out_file << "Not specified";
  }

}

/*!
 * Broadcast the AdvectDiffuse2D_InflowField variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void AdvectDiffuse2D_InflowField::Broadcast(void){
#ifdef _MPI_VERSION
  
  short i_Inflow_Field_Type_Copy;

  // Copy the value of the i_Inflow_Field_Type on the main CPU
  if (CFFC_Primary_MPI_Processor()){
    i_Inflow_Field_Type_Copy = i_Inflow_Field_Type;
  }

  // Broadcast the type of the inflow field
  MPI::COMM_WORLD.Bcast(&i_Inflow_Field_Type_Copy,
			1, 
			MPI::SHORT, 0);

  // Create the same exact solution object on each CPU different than the primary one
  // using the broadcast value.
  if (!CFFC_Primary_MPI_Processor()){
    SetInflowField(i_Inflow_Field_Type_Copy);
  }

  // Broadcast the characteristic parameters for the inflow field object
  if (IsInflowFieldSet()){
    Inflow->Broadcast();
  }
  
#endif
}

