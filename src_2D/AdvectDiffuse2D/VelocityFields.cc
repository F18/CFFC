/*!\file VelocityFields.cc
  \brief Source file initializing/implementing member variables/functions of class VelocityFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "VelocityFields.h"

short VelocityFields::i_Velocity_Field_Type = VELOCITY_FIELD_QUIESCENT; //!< still flow (no advection)
double VelocityFields::VelX = 0.0;	//!< flow velocity in x-direction set to ZERO (no advection)
double VelocityFields::VelY = 0.0;	//!< flow velocity in y-direction set to ZERO (no advection)
double VelocityFields::VelMagnitude = 0.0; //!< flow velocity magnitude set to ZERO (no advection)
double VelocityFields::theta = 0.0;	//!< flow angle set to ZERO (flow parallel to x-axis)
double VelocityFields::omega = 0.0;     //!< flow angular velocity set to ZERO (irotational flow)
Vector2D VelocityFields::CenterOfRotation = Vector2D(0.0,0.0); //!< same as the Origin
FunctionType2D VelocityFields::RotationalIntensity = VelocityFields::ConstantAngularVelocity;
short VelocityFields::i_RotationalIntensity = ANGULAR_VELOCITY_CONSTANT; //!< set to constant angular velocity


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
  } else if ( strcmp(FieldType, "Rotational_WRT_Origin") == 0 ) {
    i_Velocity_Field_Type = VELOCITY_FIELD_ROTATIONAL_WRT_ORIGIN;
  } else {
    return INVALID_INPUT_CODE;
  }
  return 0;
}

/*!
 * Set the velocity in x and y directions, based on the input parameters.
 * \param VelocityMagnitude the magnitude of the velocity vector
 * \param FlowAngle the angle of the uniform flow relative to the x-axis
 * 
 * \note the FlowAngle is in degrees!
 */
void VelocityFields::Set_UniformFlow_Parameters(const double & VelocityMagnitude,
						const double & FlowAngle){

  double rad(PI/180);

  VelX = VelocityMagnitude * cos(rad*FlowAngle);
  VelY = VelocityMagnitude * sin(rad*FlowAngle);
}

void VelocityFields::Set_RotationalFlow_Parameters(const double & ReferenceAngularVelocity,
						   FunctionType2D AngularVelocityVariation,
						   const Vector2D & ReferencePoint){

  omega = ReferenceAngularVelocity; // set the reference angular velocity
  CenterOfRotation = ReferencePoint; // set the center of rotation
  RotationalIntensity = AngularVelocityVariation; // set the function that describes the space variation of the angular velocity
}

void VelocityFields::Set_UniformRotationalFlow_Parameters(const double & AngularVelocity,
							  const Vector2D & ReferencePoint){

  Set_RotationalFlow_Parameters(AngularVelocity,ConstantAngularVelocity,ReferencePoint);
}
