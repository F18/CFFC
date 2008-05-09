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
short DiffusionFields::i_Diffusion_Field_Type_Copy = DIFFUSION_FIELD_ZERO;  //!< no diffusion
double DiffusionFields::DiffusionCoefficient = 0.0;   //!< diffusion coefficient set to ZERO (no diffusion)
double DiffusionFields::DiffusionCoefficient_Copy = 0.0;   //!< diffusion coefficient copy set to ZERO (no diffusion)
double DiffusionFields::k_x = 0.0; //!< diffusion coefficient in x-direction set to ZERO (no diffusion)
double DiffusionFields::k_y = 0.0; //!< diffusion coefficient in y-direction set to ZERO (no diffusion)
double DiffusionFields::k_u = 0.0; //!< diffusion coefficient set to ZERO (independent diffusion of the solution)
Vector2D DiffusionFields::ReferencePoint(0.0); //!< use Origin as reference point

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
  } else if ( strcmp(FieldType, "Linear") == 0 ) {
    i_Diffusion_Field_Type = DIFFUSION_FIELD_LINEAR_VARIATION;
  } else {
    return INVALID_INPUT_VALUE;
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
 * Set a linear diffusion field
 * \param  the value of the diffusion coefficient 
 */
void DiffusionFields::Set_LinearDiffusionField(const double & _k_x, const double & _k_y,
					       const double & _k_u,
					       const Vector2D & _ReferencePoint){
  i_Diffusion_Field_Type = DIFFUSION_FIELD_LINEAR_VARIATION;
  k_x = _k_x;
  k_y = _k_y;
  k_u = _k_u;
  ReferencePoint = _ReferencePoint;
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
    out_file << "\n  -> Diffusion Field: Constant, k(x,y,u) = constant"
	     << "\n     -> Diffusion Coefficient: " << DiffusionCoefficient;
    break;

  case DIFFUSION_FIELD_LINEAR_VARIATION:
    out_file << "\n  -> Diffusion Field: Linear Variation, k(x,y,u) = kx*(x-Xr) + ky*(y-Yr) + ku*u"
	     << "\n     -> kx: " << k_x
	     << "\n     -> ky: " << k_y
	     << "\n     -> ku: " << k_u
	     << "\n     -> Reference point (Xr,Yr): " << ReferencePoint;
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
  MPI::COMM_WORLD.Bcast(&k_x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&k_y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&k_u,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&ReferencePoint.y,
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
  case DIFFUSION_FIELD_LINEAR_VARIATION:
    DiffusionField = DiffusionFields::Linear_Diffusion;
    break;
  }

}
