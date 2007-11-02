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
#include "../CFD/CFD.h"

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
  static int Set_Velocity_Field_Type(const char * FieldType);
  static void Set_UniformFlow_Parameters(const double & VelocityMagnitude, const double & FlowAngle);
  static void Set_UniformRotationalFlow_Parameters(const double & AngularVelocity,
						   const Vector2D & ReferencePoint = Vector2D(0.0,0.0));
  static void Set_RotationalFlow_Parameters(const double & ReferenceAngularVelocity,
					    FunctionType2D AngularVelocityVariation = ConstantAngularVelocity,
					    const Vector2D & ReferencePoint = Vector2D(0.0,0.0));
  //@}

  //@{ @name Defined velocity fields
  static Velocity Uniform_Flow_Xdir(const double &x, const double &y);
  static Velocity Uniform_Flow_Ydir(const double &x, const double &y);
  static Velocity Uniform_Flow(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Origin(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Arbitrary_Point(const double &x, const double &y);
  //@}

  //@{ @name Defined angular velocity space variations
  /*! Imposes a constant angular velocity */
  static double ConstantAngularVelocity(const double x, const double y){ return omega; }

  /*! Imposed a space variation of the angular velocity (v_r, v_theta) = (0, omega/r), r is radial distance from center */
  static double InversVariationAngularVelocity(const double x, const double y){ return omega/sqrt(x*x + y*y); }
  //@}

  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);

  static void Print_Info(std::ostream & out_file){};

  static void Broadcast(void){};

private:
  //!< Private default constructor
  VelocityFields(void);

  //@{ @name Parameters used to determine a particular type of velocity field
  static short i_Velocity_Field_Type; //!< type indicator for the velocity field
  static double VelX;		//!< flow velocity in x-direction
  static double VelY; 		//!< flow velocity in y-direction
  static double VelMagnitude;   //!< magnitude of the velocity vector
  static double theta;		//!< flow angle (measured counterclockwise relative to the x-direction)
  static FunctionType2D RotationalIntensity;  //!< pointer to the function that describes the space variation of angular velocity
  static short i_RotationalIntensity; //!< type indicator for space variation of angular velocity
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

//! Parse the input control parameters for 
//  settings related to CENO_Execution_Mode class
template<class Input_Parameters_Type> inline
void VelocityFields::Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command){

  // Check if the next control parameter has already been identified
  if (i_command != INVALID_INPUT_CODE){
    return;
  }

  char buffer[256];

  // Try to match the next control parameter
  if (strcmp(IP.Next_Control_Parameter, "Velocity_Field") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    i_command = Set_Velocity_Field_Type(IP.Next_Control_Parameter);
  } else if (strcmp(IP.Next_Control_Parameter, "Angular_Velocity_Variation") == 0) {
    IP.Get_Next_Input_Control_Parameter();
    if ( strcmp(IP.Next_Control_Parameter, "Constant") == 0 ){
      i_RotationalIntensity = ANGULAR_VELOCITY_CONSTANT;
    } else if ( strcmp(IP.Next_Control_Parameter, "Invers_Proportional") == 0 ) {
      i_RotationalIntensity = ANGULAR_VELOCITY_INVERS_PROPORTIONAL;
    } else {
      i_command = INVALID_INPUT_CODE;
      return;
    }
    i_command = 0;
  } else if (strcmp(IP.Next_Control_Parameter, "Velocity_Xdir") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> VelX;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Velocity_Ydir") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> VelY;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Flow_Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> VelMagnitude;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Flow_Angle") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> theta;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Angular_Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> omega;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Center_Of_Rotation") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CenterOfRotation;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}


#endif
