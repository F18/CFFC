/*!\file VelocityFields.h
  \brief Header file defining 2D velocity fields. */

#ifndef _VELOCITY_FIELDS_INCLUDED
#define _VELOCITY_FIELDS_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Vector2D.h"
#include "../Utilities/TypeDefinition.h"

/*!
 * \class VelocityFields
 *
 * @brief Collection of velocity fields for the 2D problems
 *
 */

class VelocityFields{
public:
  typedef Vector2D Velocity;
  typedef Velocity (*VelocityFieldType)(const double &, const double &);

  //@{ @name Functions that set class parameters
  static void Set_UniformFlow_Velocities(const double & VelocityMagnitude, const double & FlowAngle);
  //@}

  //@{ @name Defined velocity fields
  static Velocity Uniform_Flow_Xdir(const double &x, const double &y);
  static Velocity Uniform_Flow_Ydir(const double &x, const double &y);
  static Velocity Uniform_Flow(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Origin(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Arbitrary_Point(const double &x, const double &y);
  //@}

  //@{ @name Defined angular velocity space variations
  static double ConstantAngularVelocity(const double x, const double y){ return omega; }
  static double InversVariationAngularVelocity(const double x, const double y){ return omega/sqrt(x*x + y*y); }
  //@}

private:
  //!< Private default constructor
  VelocityFields(void);

  //@{ @name Parameters used to determine a particular type of velocity field
  static double VelX;		//!< flow velocity in x-direction
  static double VelY; 		//!< flow velocity in y-direction
  static double theta;		//!< flow angle (measured counterclockwise relative to the x-direction)
  static FunctionType2D RotationalIntensity;  //!< pointer to the function that describes the space variation of angular velocity
  static double omega;		//!< flow angular velocity
  static Vector2D CenterOfRotation; //!< the center around which the flow rotates
  //@}

};

inline Vector2D VelocityFields::Uniform_Flow_Xdir(const double &x, const double &y){
  return Vector2D(VelX,0.0);
}

inline Vector2D VelocityFields::Uniform_Flow_Ydir(const double &x, const double &y){
  return Vector2D(0.0,VelY);
}

inline Vector2D VelocityFields::Uniform_Flow(const double &x, const double &y){
  return Vector2D(VelX,VelY);
}

inline Vector2D VelocityFields::Rotational_Flow_About_Origin(const double &x, const double &y){
  return RotationalIntensity(x,y)*Vector2D(-y,x);
}

inline Vector2D VelocityFields::Rotational_Flow_About_Arbitrary_Point(const double &x, const double &y){
  // do a translation of the center of rotation
  return Rotational_Flow_About_Origin(x-CenterOfRotation.x, y-CenterOfRotation.y);
}



#endif
