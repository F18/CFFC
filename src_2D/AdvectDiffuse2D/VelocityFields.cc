/*!\file VelocityFields.cc
  \brief Source file initializing/implementing member variables/functions of class VelocityFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../MPI/MPI.h"
#include "VelocityFields.h"
#include "../Utilities/Utilities.h"

short  VelocityFields::i_Velocity_Field_Type = VELOCITY_FIELD_QUIESCENT;     //!< still flow (no advection)
VelocityFields::Velocity VelocityFields::CartesianVelocity = Vector2D(0.0);  //!< flow velocity set to ZERO (no advection)
Vector2D VelocityFields::MagnitudeAngleVelocity = Vector2D(0.0);  //!< flow velocity set to ZERO (no advection)
double VelocityFields::omega = 0.0;                               //!< angular flow velocity set to ZERO (irotational flow)
Vector2D VelocityFields::CenterOfRotation = Vector2D(0.0,0.0);    //!< same as the Origin
FunctionType2D VelocityFields::RotationalIntensity = VelocityFields::ConstantAngularVelocity;
short VelocityFields::i_RotationalIntensity = ANGULAR_VELOCITY_CONSTANT;  //!< set to constant angular velocity


// ===  Member functions ===
/*!
 * Set the velocity field type index based on the input 
 * \param FieldType C-type string providing the name of the velocity field
 * \return If the FieldType is recognized the returned value is '0', otherwise is 'INVALID_INPUT_CODE'
 */
int VelocityFields::Set_Velocity_Field_Type(const char * FieldType){
  if ( strcmp(FieldType, "Quiescent") == 0 ){
    i_Velocity_Field_Type = VELOCITY_FIELD_QUIESCENT;
  } else if ( strcmp(FieldType, "Uniform") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM;
  } else if ( strcmp(FieldType, "Uniform_Xdir") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_X_DIRECTION;
  } else if ( strcmp(FieldType, "Uniform_Ydir") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_Y_DIRECTION;
  } else if ( strcmp(FieldType, "Rotational") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_ROTATIONAL_WRT_ARBITRATY_POINT;
  } else if ( strcmp(FieldType, "Rotational_Uniform") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_ROTATIONAL_WRT_ARBITRATY_POINT;
  } else {
    return INVALID_INPUT_CODE;
  }
  return 0;
}

/*!
 * Set the corresponding spatial variation of the angular velocity
 * \param VariationType C-type string providing the name of the angular velocity variation
 * \return If the VariationType is recognized the returned value is '0', otherwise is 'INVALID_INPUT_CODE'
 */
int VelocityFields::Set_Angular_Velocity_Variation(const char * VariationType){
  if ( strcmp(VariationType, "Constant") == 0 ){
    RotationalIntensity = ConstantAngularVelocity;
    i_RotationalIntensity = ANGULAR_VELOCITY_CONSTANT;
  } else if ( strcmp(VariationType, "Inverse_Proportional_Distance") == 0 ) {
    RotationalIntensity = InverseVariationAngularVelocity;
    i_RotationalIntensity = ANGULAR_VELOCITY_INVERSE_PROPORTIONAL;
  } else {
    return INVALID_INPUT_CODE;
  }
  return 0;
}

/*!
 * Set the RotationalIntensity pointer to the spatial variation of the angular velocity
 * \param VariationIndex the index of the angular velocity variation
 * \throw runtime_error if the index is not recognized
 */
void VelocityFields::Set_Angular_Velocity_Variation(void){
  switch(i_RotationalIntensity){
  case ANGULAR_VELOCITY_CONSTANT:
    RotationalIntensity = ConstantAngularVelocity;
    break;
    
  case ANGULAR_VELOCITY_INVERSE_PROPORTIONAL:
    RotationalIntensity = InverseVariationAngularVelocity;
    break;

  default:
    throw runtime_error("VelocityFields::Set_Angular_Velocity_Variation() ERROR: Unknown angular velocity variation index");
  }
}

/*!
 * Compute Cartesian velocity components based on the 
 * magnitude of the velocity vector and the flow angle
 */
void VelocityFields::ComputeCartesianVelocity(void){
  // x-axis component
  CartesianVelocity.x = MagnitudeAngleVelocity.x * cos(MagnitudeAngleVelocity.y);
  // y-axis component
  CartesianVelocity.y = MagnitudeAngleVelocity.x * sin(MagnitudeAngleVelocity.y);
}

/*!
 * Compute magnitude of the velocity vector and the 
 * flow angle based on the Cartesian velocity.
 */
void VelocityFields::ComputeMagnitudeAngleVelocity(void){
  // magnitude of the velocity vector
  MagnitudeAngleVelocity.x = sqrt(CartesianVelocity.x*CartesianVelocity.x + CartesianVelocity.y*CartesianVelocity.y);
  // flow angle
  MagnitudeAngleVelocity.y = arctan(CartesianVelocity.y,CartesianVelocity.x);
}

/*!
 * Set a uniform flow which has the flow velocity provided by the input parameters.
 * \param VelocityMagnitude the magnitude of the velocity vector
 * \param FlowAngle the angle of the uniform flow relative to the x-axis
 * 
 * \note FlowAngle is in degrees!
 */
void VelocityFields::Set_UniformFlow_MagnitudeAngleVelocity(const double & VelocityMagnitude,
							    const double & FlowAngle){

  // set the flow type
  if (VelocityMagnitude == 0.0){
    i_Velocity_Field_Type = VELOCITY_FIELD_QUIESCENT;
  } else if ( FlowAngle == 0.0 || FlowAngle == 180){
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_X_DIRECTION;
  } else if (FlowAngle == 90.0 || FlowAngle == 270){
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_Y_DIRECTION;
  } else {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM;
  }

  // set the MagnitudeAngleVelocity 
  MagnitudeAngleVelocity.x = VelocityMagnitude;
  MagnitudeAngleVelocity.y = Radians(FlowAngle);

  // set the CartesianVelocity
  ComputeCartesianVelocity();
}

/*!
 * Set a uniform flow which has the flow velocity provided by the input parameters.
 * \param Velocity_XAxis velocity in the x-direction
 * \param Velocity_YAxis velocity in the y-direction
 */
void VelocityFields::Set_UniformFlow_CartesianVelocity(const double & Velocity_XAxis,
						       const double & Velocity_YAxis){

  // set the flow type
  if (Velocity_XAxis == 0.0 && Velocity_YAxis == 0.0){
    i_Velocity_Field_Type = VELOCITY_FIELD_QUIESCENT;
  } else if (Velocity_XAxis != 0.0 && Velocity_YAxis == 0.0){
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_X_DIRECTION;
  } else if (Velocity_XAxis == 0.0 && Velocity_YAxis != 0.0){
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_Y_DIRECTION;
  } else {
    i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM;
  }

  // set the CartesianVelocity 
  CartesianVelocity.x = Velocity_XAxis;
  CartesianVelocity.y = Velocity_YAxis;

  // set the MagnitudeAngleVelocity 
  ComputeMagnitudeAngleVelocity();
}

/*!
 * Set a rotational flow with a variational angular flow velocity
 * \param ReferenceAngularVelocity the reference angular flow velocity
 * \param AngularVelocityVariation C-type string providing the name of the angular velocity variation
 * \param ReferencePoint the center of rotation
 */
void VelocityFields::Set_RotationalFlow_Parameters(const double & ReferenceAngularVelocity,
						   const char * AngularVelocityVariation,
						   const Vector2D & ReferencePoint){

  // set the flow type
  i_Velocity_Field_Type = VELOCITY_FIELD_ROTATIONAL_WRT_ARBITRATY_POINT;

  omega = ReferenceAngularVelocity;  // set the reference angular velocity
  CenterOfRotation = ReferencePoint; // set the center of rotation

  // set the function that describes the space variation of the angular velocity
  Set_Angular_Velocity_Variation(AngularVelocityVariation); 
}

/*!
 * Set a rotational flow with a constant angular flow velocity
 * \param AngularVelocity the angular flow velocity
 * \param ReferencePoint the center of rotation
 */
void VelocityFields::Set_UniformRotationalFlow_Parameters(const double & AngularVelocity,
							  const Vector2D & ReferencePoint){

  Set_RotationalFlow_Parameters(AngularVelocity,"Constant",ReferencePoint);

  // reset the flow type
  i_Velocity_Field_Type = VELOCITY_FIELD_UNIFORM_ROTATIONAL_WRT_ARBITRATY_POINT;
}

/*!
 * Print the relevant parameters of the VelocityFields class for the 
 * type of selected velocity field to the provided output stream.
 */
void VelocityFields::Print_Info(std::ostream & out_file){

  switch(i_Velocity_Field_Type){
  case VELOCITY_FIELD_QUIESCENT:
    out_file << "\n  -> Flow Velocity Field: Quiescent Flow";
    break;

  case VELOCITY_FIELD_UNIFORM:
    out_file << "\n  -> Flow Velocity Field: Uniform Flow";
    out_file << "\n     -> Velocity Magnitude: " << MagnitudeAngleVelocity.x;
    out_file << "\n     -> Flow Angle: " << Degrees(MagnitudeAngleVelocity.y) << " degrees";
    break;

  case VELOCITY_FIELD_UNIFORM_X_DIRECTION:
    out_file << "\n  -> Flow Velocity Field: Uniform Flow";
    out_file << "\n     -> Velocity Magnitude: " << CartesianVelocity.x;
    if (CartesianVelocity.x > 0){
      out_file << "\n     -> Flow Angle: " << 0 << " degrees (parallel to x-axis)";
    } else if (CartesianVelocity.x < 0){
      out_file << "\n     -> Flow Angle: " << 180 << " degrees (parallel to x-axis)";
    }
    break;

  case VELOCITY_FIELD_UNIFORM_Y_DIRECTION:
    out_file << "\n  -> Flow Velocity Field: Uniform Flow";
    out_file << "\n     -> Velocity Magnitude: " << CartesianVelocity.y;
    if (CartesianVelocity.y > 0){
      out_file << "\n     -> Flow Angle: " << 90 << " degrees (parallel to y-axis)";
    } else if (CartesianVelocity.y < 0){
      out_file << "\n     -> Flow Angle: " << 270  << " degrees (parallel to y-axis)";
    }
    break;

  case VELOCITY_FIELD_ROTATIONAL_WRT_ARBITRATY_POINT:
    out_file << "\n  -> Flow Velocity Field: Rotationary Flow";
    out_file << "\n     -> Center of Rotation: " << CenterOfRotation;
    out_file << "\n     -> Angular Velocity Variation: ";
    switch(i_RotationalIntensity){
    case ANGULAR_VELOCITY_CONSTANT:
      out_file << "Constant";
      out_file << "\n     -> Angular velocity: " << omega;
      break;
    case ANGULAR_VELOCITY_INVERSE_PROPORTIONAL:
      out_file << "Inverse-proportional distance";
      out_file << "\n     -> Reference angular velocity: " << omega;
      break;
    } //endswitch
    break;

  case VELOCITY_FIELD_UNIFORM_ROTATIONAL_WRT_ARBITRATY_POINT:
    out_file << "\n  -> Flow Velocity Field: Uniform Rotationary Flow";
    out_file << "\n     -> Center of Rotation: " << CenterOfRotation;
    // make sure that the RotationalIntensity pointer is set correctly
    if (i_RotationalIntensity != ANGULAR_VELOCITY_CONSTANT){
      RotationalIntensity = ConstantAngularVelocity;
      i_RotationalIntensity = ANGULAR_VELOCITY_CONSTANT;
    }
    out_file << "\n     -> Angular Velocity Variation: Constant";
    out_file << "\n     -> Angular velocity: " << omega;
    break;
  }  
}

/*!
 * Broadcast the VelocityFields variables to all      
 * processors associated with the specified communicator
 * from the specified processor using the MPI broadcast 
 * routine.
 *
 * \todo Switch to a user-difined datatype
 */
void VelocityFields::Broadcast(void){
#ifdef _MPI_VERSION
  
  MPI::COMM_WORLD.Bcast(&i_Velocity_Field_Type,
			1, 
			MPI::SHORT, 0);
  MPI::COMM_WORLD.Bcast(&CartesianVelocity.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&CartesianVelocity.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&MagnitudeAngleVelocity.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&MagnitudeAngleVelocity.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&omega,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&CenterOfRotation.x,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&CenterOfRotation.y,
			1, 
			MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&i_RotationalIntensity,
			1, 
			MPI::SHORT, 0);
  // Set RotationalIntensity pointer based on i_RotationalIntensity value
  Set_Angular_Velocity_Variation();
  
#endif
}

/*!
 * Set the passed pointer to the current velocity flow field.
 *
 * \param VelocityField this pointer will be set to point to the current flow field
 */
void VelocityFields::Connect_Pointer_To_Flow_Field(VelocityFields::VelocityFieldTypeRef VelocityField){

  switch(i_Velocity_Field_Type){
  case VELOCITY_FIELD_QUIESCENT:
    VelocityField = VelocityFields::Quiescent_Flow;
    break;
  case VELOCITY_FIELD_UNIFORM:
    VelocityField = VelocityFields::Uniform_Flow;
    break;
  case VELOCITY_FIELD_UNIFORM_X_DIRECTION:
    VelocityField = VelocityFields::Uniform_Flow_XAxis;
    break;
  case VELOCITY_FIELD_UNIFORM_Y_DIRECTION:
    VelocityField = VelocityFields::Uniform_Flow_YAxis;
    break;
  case VELOCITY_FIELD_ROTATIONAL_WRT_ARBITRATY_POINT:
    VelocityField = VelocityFields::Rotational_Flow_About_Arbitrary_Point;
    break;
  case VELOCITY_FIELD_UNIFORM_ROTATIONAL_WRT_ARBITRATY_POINT:
    VelocityField = VelocityFields::Rotational_Flow_About_Arbitrary_Point;
    break;
  }

}
