/*!\file SourceTermFields.cc
  \brief Source file initializing/implementing member variables/functions of class SourceTermFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "SourceTermFields.h"
#include "AdvectDiffuse2DInput.h"

// ==== Member variables ====
SourceFieldBasicType * SourceTermFields::SourcePtr = NULL; //!< no associated source field
short SourceTermFields::i_Source_Field_Type = -1;          //!< value for no associated source field

// ==== Member functions ====
/*! Default constructor is used to create the unique object of this class */
SourceTermFields::SourceTermFields(void){ }

/*! Destroy the source field object */
void SourceTermFields::DestroySource(void){
 delete SourcePtr;
 SourcePtr = NULL;		// set to no associated source field
 i_Source_Field_Type = -1;	// set to value for no associated source field
}

/*!
 * Create and return the address of the source term field.
 * This is a singleton class (only one instance is created).\n
 * DON'T modify this function unless you are familiar with the Meyers singleton pattern.
 */
SourceTermFields& SourceTermFields::getInstance(void){
  static SourceTermFields inst;
  std::atexit(DestroySource);		// schedule for deallocation at exit
  return inst;
}

/*!
 * Set the source field based on the required field type
 */
void SourceTermFields::SetSourceField(const short &FieldIndex){
  // deallocate the source term
  DestroySource();

  // create the proper source field and set the index accordingly
  switch (FieldIndex){
  case SOURCE_FIELD_ZERO:
    SourcePtr = new ZERO_SourceField;
    break;
    
  case SOURCE_FIELD_LINEAR_VARIATION:
    SourcePtr = new Linear_SourceField;
    break;
    
  case SOURCE_FIELD_EXPONENTIAL_VARIATION:
    SourcePtr = new Exponential_SourceField;
    break;

  default:
    throw runtime_error("SourceTermFields::SetSourceField() ERROR! Unknown source term field type.");
  }

  // Store the field type
  i_Source_Field_Type = FieldIndex;
}

/*!
 * Parse the input control parameters for settings 
 * related to SourceTermFields class.
 * The first entry related to the source term in the
 * input file MUST specify the type of the source term.
 * Consequent parameters will be parsed with the parser 
 * particular to the specified source term type.
 */
void SourceTermFields::Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP,
							  int & i_command){
    // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Source_Field") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "No_Source") == 0 ){
      SetSourceField(SOURCE_FIELD_ZERO);
    } else if ( strcmp(IP.Next_Control_Parameter, "Linear_Variation") == 0 ) {
      SetSourceField(SOURCE_FIELD_LINEAR_VARIATION);
    } else if ( strcmp(IP.Next_Control_Parameter, "Exponential_Variation") == 0 ) {
      SetSourceField(SOURCE_FIELD_EXPONENTIAL_VARIATION);
    } else {
      i_command = INVALID_INPUT_VALUE;
      return;
    }
    i_command = 0;

  } else {
    // Continue parsing with the parser of the current source field
    if (SourcePtr != NULL){
      SourcePtr->Parse_Next_Input_Control_Parameter(IP,i_command);
    }
    return;
  } // endif
}

/*!
 * Print the relevant parameters of the SourceTermFields class for the 
 * selected source field type to the provided output stream.
 */
void SourceTermFields::Print_Info(std::ostream & out_file){

  out_file << "\n  -> Source Field: " ; 
  if (SourcePtr != NULL){
    // Output field name
    out_file << SourcePtr->whatField();

    // Output field characteristic parameters
    SourcePtr->Print_Info(out_file);
  } else {
    // There is no associated source field
    out_file << "Not specified";
  }

}

/*!
 * Broadcast the SourceTermFields variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 */
void SourceTermFields::Broadcast(void){
#ifdef _MPI_VERSION
  
  short i_Source_Field_Type_Copy;

  // Copy the value of the i_Source_Field_Type on the main CPU
  if (CFFC_Primary_MPI_Processor()){
    i_Source_Field_Type_Copy = i_Source_Field_Type;
  }

  // Broadcast the type of the source field
  MPI::COMM_WORLD.Bcast(&i_Source_Field_Type_Copy,
			1, 
			MPI::SHORT, 0);

  // Create the same source field object on each CPU different than the primary one
  // using the broadcast value.
  if (!CFFC_Primary_MPI_Processor()){
    SetSourceField(i_Source_Field_Type_Copy);
  }

  // Broadcast the characteristic parameters for the source field object
  SourcePtr->Broadcast();
  
#endif
}

/*!
 * Determine the stability limit imposed by the 
 * non-linear source term
 */
double SourceTermFields::getStabilityLimit(const double &x, const double &y, const double &u){
  // Pass possibly necessary information for determining the stability limit

  // Return the absolute value to make sure that DeltaT gets set correctly
  if (FieldRequireIntegration()) {
    return fabs(SourcePtr->StabilityLimit(x,y,u));
  } else {
    return fabs(SourcePtr->StabilityLimit(u));
  }
}
