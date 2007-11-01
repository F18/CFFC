/*!\file VelocityFields.cc
  \brief Source file initializing/implementing member variables/functions of class VelocityFields. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "VelocityFields.h"

double VelocityFields::VelX = 0.0;	//!< flow velocity in x-direction set to ZERO (no advection)
double VelocityFields::VelY = 0.0;	//!< flow velocity in y-direction set to ZERO (no advection)
double VelocityFields::theta = 0.0;	//!< flow angle set to ZERO (flow parallel to x-axis)
double VelocityFields::omega = 0.0;     //!< flow angular velocity set to ZERO (irotational flow)
Vector2D VelocityFields::CenterOfRotation = Vector2D(0.0,0.0); //!< same as the Origin
FunctionType2D VelocityFields::RotationalIntensity = VelocityFields::ConstantAngularVelocity;

// ===  Member functions ===
void VelocityFields::Set_UniformFlow_Velocities(const double & VelocityMagnitude,
						const double & FlowAngle){
  VelX = VelocityMagnitude * cos(FlowAngle);
  VelY = VelocityMagnitude * sin(FlowAngle);
}
