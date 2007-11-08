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
#include "../CFD/CFD.h"

/*!
 * \class VelocityFields
 *
 * @brief Collection of velocity fields for 2D problems
 * \nosubgrouping
 */

class VelocityFields{
public:
  //! @name Defined public types
  //@{
  typedef Vector2D Velocity;	//!< Velocity type
  typedef Velocity (*VelocityFieldType)(const double &, const double &); //!< Velocity field type
  typedef Velocity (*&VelocityFieldTypeRef)(const double &, const double &); //!< Reference to velocity field type
  //@}

  //! @name Functions that set class parameters
  //@{
  static void Set_UniformFlow_MagnitudeAngleVelocity(const double & VelocityMagnitude, const double & FlowAngle);
  static void Set_UniformFlow_CartesianVelocity(const double & Velocity_XAxis = 0.0, const double & Velocity_YAxis = 0.0);
  static void Set_UniformRotationalFlow_Parameters(const double & AngularVelocity,
						   const Vector2D & ReferencePoint = Vector2D(0.0,0.0));
  static void Set_RotationalFlow_Parameters(const double & ReferenceAngularVelocity,
					    const char * AngularVelocityVariation = "Constant",
					    const Vector2D & ReferencePoint = Vector2D(0.0,0.0));
  //@}

  //! @name Defined velocity fields
  //@{
  static Velocity Quiescent_Flow(const double &x, const double &y);
  static Velocity Uniform_Flow_XAxis(const double &x, const double &y);
  static Velocity Uniform_Flow_YAxis(const double &x, const double &y);
  static Velocity Uniform_Flow(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Origin(const double &x, const double &y);
  static Velocity Rotational_Flow_About_Arbitrary_Point(const double &x, const double &y);
  //@}

  //! @name Defined angular velocity space variations
  //@{
  /*! Imposes a constant angular velocity */
  static double ConstantAngularVelocity(const double x, const double y){ return omega; }

  /*! Imposed a space variation of the angular velocity (v_r, v_theta) = (0, omega/r), r is radial distance from center */
  static double InverseVariationAngularVelocity(const double x, const double y){ return omega/sqrt(x*x + y*y); }
  //@}

  //! @name Functions for input-output and broadcast
  //@{
  template<class Input_Parameters_Type>
  static void Parse_Next_Input_Control_Parameter(Input_Parameters_Type & IP, int & i_command);
  static void Print_Info(std::ostream & out_file);
  static void Broadcast(void);
  //@}

  //! @name Functions for setting external pointers
  //@{
  static void Connect_Pointer_To_Flow_Field(VelocityFieldTypeRef VelocityField);
  //@}

protected:
  VelocityFields(void);   //!< Private default constructor
  VelocityFields(const VelocityFields&); //!< Private copy constructor
  VelocityFields& operator=(const VelocityFields&); //!< Private assignment operator

  //! @name Parameters used to determine a particular type of velocity field
  //@{
  static short i_Velocity_Field_Type;     //!< type indicator for the velocity field
  static Velocity CartesianVelocity;      //!< flow velocity in the x-axis and y-axis
  static Vector2D MagnitudeAngleVelocity; //!< velocity given by magnitude and angle relative 
                                          // to the x-axis (measured counterclockwise in radians)
  static FunctionType2D RotationalIntensity;  //!< pointer to the function that describes the space variation of angular velocity
  static short i_RotationalIntensity;         //!< type indicator for space variation of angular velocity
  static double omega;		              //!< flow angular velocity
  static Vector2D CenterOfRotation;           //!< the center around which the flow rotates
  //@}

  //! @name Transformation functions
  //@{
  static void ComputeCartesianVelocity(void);
  static void ComputeMagnitudeAngleVelocity(void);
  static double Radians(const double &AngleInDegrees){return (PI*AngleInDegrees)/180.0; }
  static double Degrees(const double &AngleInRadians){return (180*AngleInRadians)/PI; }
  //@}

  //! @name Set field types
  //@{
  static int Set_Velocity_Field_Type(const char * FieldType);
  static int Set_Angular_Velocity_Variation(const char * VariationType);
  static void Set_Angular_Velocity_Variation(void);
  //@}

};

/*! Defines a quiescent flow (zero velocity) */
inline Vector2D VelocityFields::Quiescent_Flow(const double &x, const double &y){
  return Vector2D(0.0);
}

/*! Defines a uniform flow parallel with x-axis */
inline Vector2D VelocityFields::Uniform_Flow_XAxis(const double &x, const double &y){
  return Vector2D(CartesianVelocity.x,0.0);
}

/*! Defines a uniform flow parallel with y-axis */
inline Vector2D VelocityFields::Uniform_Flow_YAxis(const double &x, const double &y){
  return Vector2D(0.0,CartesianVelocity.y);
}

/*! Defines a uniform flow */
inline Vector2D VelocityFields::Uniform_Flow(const double &x, const double &y){
  return Vector2D(CartesianVelocity.x,CartesianVelocity.y);
}

/*!
 * Defines a rotational flow centered in (0,0) (Origin) that has an angular velocity given by the RotationalIntensity function
 */
inline Vector2D VelocityFields::Rotational_Flow_About_Origin(const double &x, const double &y){
  return RotationalIntensity(x,y)*Vector2D(-y,x);
}

/*!
 * Defines a rotational flow centered in an arbitrary point that has an angular velocity given by the RotationalIntensity function
 */
inline Vector2D VelocityFields::Rotational_Flow_About_Arbitrary_Point(const double &x, const double &y){
  // do a translation of the center of rotation
  return Rotational_Flow_About_Origin(x-CenterOfRotation.x, y-CenterOfRotation.y);
}

//! Parse the input control parameters for 
//  settings related to VelocityFields class
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
    i_command = Set_Angular_Velocity_Variation(IP.Next_Control_Parameter);
  } else if (strcmp(IP.Next_Control_Parameter, "Velocity_Xdir") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CartesianVelocity.x;
    IP.Input_File.getline(buffer, sizeof(buffer));
    ComputeMagnitudeAngleVelocity();	// update the polar velocities
  } else if (strcmp(IP.Next_Control_Parameter, "Velocity_Ydir") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CartesianVelocity.y;
    IP.Input_File.getline(buffer, sizeof(buffer));
    ComputeMagnitudeAngleVelocity();	// update the polar velocities
  } else if (strcmp(IP.Next_Control_Parameter, "Flow_Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> MagnitudeAngleVelocity.x;
    IP.Input_File.getline(buffer, sizeof(buffer));
    ComputeCartesianVelocity();	// update the polar velocities
  } else if (strcmp(IP.Next_Control_Parameter, "Flow_Angle") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> MagnitudeAngleVelocity.y;
    // transform the degrees to radians
    MagnitudeAngleVelocity.y = Radians(MagnitudeAngleVelocity.y);
    IP.Input_File.getline(buffer, sizeof(buffer));
    ComputeCartesianVelocity();	// update the polar velocities
  } else if (strcmp(IP.Next_Control_Parameter, "Angular_Velocity") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> omega;
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else if (strcmp(IP.Next_Control_Parameter, "Center_Of_Rotation") == 0) {
    i_command = 0;
    ++IP.Line_Number;
    IP.Input_File >> CenterOfRotation;
    IP.Input_File.setf(ios::skipws);
    IP.Input_File.getline(buffer, sizeof(buffer));
  } else {
    i_command = INVALID_INPUT_CODE;
  } // endif

}


#endif
