/*!\file DiffusionFields.cc
  \brief Source file initializing/implementing member variables/functions of class DiffusionFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "DiffusionFields.h"
#include "../Utilities/Utilities.h"

// ===  Static member variables ===
short DiffusionFields::i_Diffusion_Field_Type = DIFFUSION_FIELD_ZERO;  //!< no diffusion
double DiffusionFields::DiffusionCoefficient = 0.0;   //!< diffusion coefficient set to ZERO (no diffusion)

// ===  Member functions ===
/*!
 * Set the non-linear diffusion field type index based on the input 
 * \param FieldType C-type string providing the name of the diffusion field
 * \return If the FieldType is recognized the returned value is '0', otherwise is 'INVALID_INPUT_CODE'
 */
int DiffusionFields::Set_Diffusion_Field_Type(const char * FieldType){
  if ( strcmp(FieldType, "No_Diffusion") == 0 ){
    i_Diffusion_Field_Type = DIFFUSION_FIELD_ZERO;
  } else if ( strcmp(FieldType, "Constant") == 0 ) {
    i_Diffusion_Field_Type = DIFFUSION_FIELD_CONSTANT;
  } else {
    return INVALID_INPUT_CODE;
  }
  return 0;
}

/*!
 * Set a constant diffusion field
 * \param DiffusionCoeff the value of the diffusion coefficient 
 */
void DiffusionFields::Set_ConstantDiffusionField(const double & DiffusionCoeff){
  i_Diffusion_Field_Type = DIFFUSION_FIELD_CONSTANT;
  DiffusionCoefficient = DiffusionCoeff;
}

/*!
 * Print the relevant parameters of the DiffusionFields class for the 
 * type of selected non-linear diffusion field to the provided output stream.
 */
void DiffusionFields::Print_Info(std::ostream & out_file){

  switch(i_Diffusion_Field_Type){
  case DIFFUSION_FIELD_ZERO:
    out_file << "\n  -> Diffusion Field: No Diffusion";
    break;

  case DIFFUSION_FIELD_CONSTANT:
    out_file << "\n  -> Diffusion Field: Constant";
    out_file << "\n     -> Diffusion Coefficient: " << DiffusionCoefficient;
    break;
  }
}

/*!
 * Broadcast the DiffusionFields variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 *
 * \todo Switch to a user-defined datatype
 */
void DiffusionFields::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&i_Diffusion_Field_Type,
			1, 
			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&DiffusionCoefficient,
			1, 
			MPI::DOUBLE, 0);
  
#endif
}

/*!
 * Set the passed pointer to the current non-linear diffusion field.
 *
 * \param DiffusionField this pointer will be set to point to the current diffusion field
 */
void DiffusionFields::Connect_Pointer_To_Diffusion_Field(DiffusionFields::NonlinearDiffusionFieldTypeRef DiffusionField){

  switch(i_Diffusion_Field_Type){
  case DIFFUSION_FIELD_ZERO:
    DiffusionField = DiffusionFields::Zero_Diffusion;
    break;
  case DIFFUSION_FIELD_CONSTANT:
    DiffusionField = DiffusionFields::Constant_Diffusion;
    break;
  }

}

