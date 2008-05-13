/*!\file InflowFields.h
  \brief Header file defining 2D inflow fields. */

#ifndef _INFLOWFIELDS_INCLUDED
#define _INFLOWFIELDS_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Math/Math.h"
#include "../Math/Vector2D.h"
#include "../Math/NumericalLibrary.h"
#include "../Utilities/TypeDefinition.h"
#include "../Utilities/Utilities.h"
#include "../CFD/CFD.h"


/******************************************
 *                                        *
 *         Inflow Field Types             *
 *                                        *
 *****************************************/
#define NO_INFLOW_FIELD             -1
#define SINUSOIDAL_I                 0
#define SINUSOIDAL_II                1
#define SINUSOIDAL_III               2
#define SINUSOIDAL_IV                3
#define CONSTANT_INFLOW_FIELD        4
#define HYPERBOLIC_TANGENT_I         5
#define EXPONENTIAL_SINUSOIDAL       6



// Declare the input parameters class
class AdvectDiffuse2D_Input_Parameters;

/*! 
 * \class InflowFieldBasicType
 *
 * \brief Basic data type for any inflow field.
 *
 * This is an abstract data type (ADT).
 */
class InflowFieldBasicType{
public:

  //! Default ctor
  InflowFieldBasicType(void);

  //! @name Virtual member functions
  //@{
  //! Declare a pure virtual destructor
  virtual ~InflowFieldBasicType(void) = 0;

  /*! Calculate the inflow field value at the location of interest */
  virtual double EvaluateSolutionAt(const double &x, const double &y) const = 0;

  //! Update internal variables
  virtual void Set_InflowField_Parameters(void){ };

  //! @name Functions for input-output and broadcast
  //@{
  virtual void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command) = 0;
  virtual void Print_Info(std::ostream & out_file) = 0;
  virtual void Broadcast(void) = 0;
  //@}
  //@}

  const string & whatInflowField(void) const { return InflowFieldName; }   //!< Get the name of the inflow field

protected:
  string InflowFieldName;		//!< storage for the name of the inflow field

};


/*! 
 * \class Sinusoidal_I_InflowField
 * 
 * \brief Implements a sinusoidal inflow field with the following expression: 
 *        \f$ Inflow(x,y) = \sin(\pi \, y) \f$
 */
class Sinusoidal_I_InflowField: public InflowFieldBasicType{
public:

  //! Basic Constructor
  Sinusoidal_I_InflowField(void){ 
    InflowFieldName = "Sinusoidal 1, Inflow(x,y) = sin(PI*y)";	// Name the inflow field
  }

  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const {return sin(PI*y);}

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command){ };

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file){ };

  //! Broadcast relevant parameters
  void Broadcast(void){ };

private:
};


/*! 
 * \class Sinusoidal_II_InflowField
 * 
 * \brief Implements a sinusoidal inflow field with the following expression: 
 *        \f$ Inflow(x,y) = \sin \left(\frac{\pi \ln r}{\ln 2} \right) \f$ \n
 *
 * where, r is the distance between the location of interest and the reference point
 */
class Sinusoidal_II_InflowField: public InflowFieldBasicType{
public:

  //! Basic Constructor
  Sinusoidal_II_InflowField(void): ReferencePoint(0.0){ 
    InflowFieldName = "Sinusoidal 2, Inflow(x,y) = sin(PI*ln(r)/ln(2))";	// Name the inflow field
  }

  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const {
    return sin(PI*log(abs(Vector2D(x,y)-ReferencePoint))/log(2));
  }

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);

private:
  Vector2D ReferencePoint;	//!< 2D-space location used as reference for the field
};


/*! 
 * \class Sinusoidal_III_InflowField
 * 
 * \brief Implements a sinusoidal inflow field with the following expression: 
 *        \f$ Inflow(x,y) = \sin \left(\pi \frac{r - A}{B} \right) \f$ \n
 *
 * if r is between [R_Min,R_Max], where r is the distance between the location 
 * of interest and the reference point. \n
 * For r outside of the domain, the inflow is 0.0.
 */
class Sinusoidal_III_InflowField: public InflowFieldBasicType{
public:
  
  //! Basic Constructor
  Sinusoidal_III_InflowField(void): ReferencePoint(0.0),
				    A(0.1), B(0.6),
				    R_Min(0.0), R_Max(1.0){ 

    InflowFieldName = "Sinusoidal 3, Inflow(x,y) = sin(PI*(r-A)/B), between [R_Min,R_Max]";	// Name the inflow field
  }
  
  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);
  
private:
  Vector2D ReferencePoint;	//!< 2D-space location used as reference for the field
  double A,B;			//!< Constant parameters of the inflow field
  double R_Min, R_Max;		//!< Minimum and maximum values for r
};

//! Return value of the inflow field
inline double Sinusoidal_III_InflowField::EvaluateSolutionAt(const double &x, const double &y) const {

  double R;

  // Calculate the distance R relative to the reference point
  R = abs(Vector2D(x,y) - ReferencePoint);

  if (R < R_Min || R > R_Max){
    return 0.0;
  } else {
    return sin(PI*(R-A)/B);
  }
}


/*! 
 * \class Sinusoidal_IV_InflowField
 * 
 * \brief Implements a sinusoidal inflow field with the following expression: 
 *        \f$ Inflow(x,y) = \sin^2 \left(\pi \frac{r - A}{1 - 2A} \right) \f$ \n
 *
 * if r is between [R_Min,R_Max], where r is the distance between the location 
 * of interest and the reference point. \n
 * For r outside of the domain, the inflow is 0.0.
 */
class Sinusoidal_IV_InflowField: public InflowFieldBasicType{
public:
  
  //! Basic Constructor
  Sinusoidal_IV_InflowField(void): ReferencePoint(0.0),
				   A(0.2),
				   R_Min(0.2), R_Max(0.8){ 
    InflowFieldName = "Sinusoidal 4, Inflow(x,y) = sin^2(PI*(r-A)/(1-2A)), between [R_Min,R_Max]";  // Name the inflow field
  }
  
  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);
  
private:
  Vector2D ReferencePoint;	//!< 2D-space location used as reference for the field
  double A;			//!< Constant parameter of the inflow field
  double R_Min, R_Max;		//!< Minimum and maximum values for r
};

//! Return value of the inflow field
inline double Sinusoidal_IV_InflowField::EvaluateSolutionAt(const double &x, const double &y) const {

  double R;

  // Calculate the distance R relative to the reference point
  R = abs(Vector2D(x,y) - ReferencePoint);

  if (R < R_Min || R > R_Max){
    return 0.0;
  } else {
    return sqr(sin(PI*(R-A)/(1.0-2.0*A)));
  }
}


/*! 
 * \class Constant_InflowField
 * 
 * \brief Implements a constant value inflow field with the following expression: 
 *        \f$ Inflow(x,y) = A \f$ \n
 */
class Constant_InflowField: public InflowFieldBasicType{
public:
  
  //! Basic Constructor
  Constant_InflowField(void): A(1.0){ 
    InflowFieldName = "Constant, Inflow(x,y) = A";  // Name the inflow field
  }
  
  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const { return A; };

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);
  
private:
  double A;			//!< Constant parameter of the inflow field
};

/*! 
 * \class Hyperbolic_Tangent_I_InflowField
 * 
 * \brief Implements an inflow field with the following expression: 
 *        \f$ Inflow(x,y) = M \left\{ A - \tanh^2 [S (r - r_{0}) ] \right\} \f$ \n
 *
 * where r is the location of interest, \f$r_{0}\f$ is the reference point,
 * M is the magnitude of the function, A is a shift coefficient and S controls the
 * steepness. \n
 */
class Hyperbolic_Tangent_I_InflowField: public InflowFieldBasicType{
public:
  
  //! Basic Constructor
  Hyperbolic_Tangent_I_InflowField(void): ReferencePoint(0.0),
					  A(1.0), Magnitude(1.0),
					  Steepness(1.0)
  { 

    InflowFieldName = "Hyperbolic Tan 1, Inflow(x,y) = M*[A - tanh^2[S*(r - r0)] ]";	// Name the inflow field
  }
  
  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);
  
private:
  Vector2D ReferencePoint;	//!< 2D-space location used as reference for the field
  double A;			//!< Constant parameter of the inflow field
  double Magnitude, Steepness; 	//!< Function magnitude and steepness control
};

//! Return value of the inflow field
inline double Hyperbolic_Tangent_I_InflowField::EvaluateSolutionAt(const double &x, const double &y) const {

  double R;
  double Func;

  // Calculate the distance R relative to the reference point
  R = abs(Vector2D(x,y) - ReferencePoint);

  Func  = tanh(Steepness * R);
  Func *= Func;

  return Magnitude * (A - Func); 
}

/*! 
 * \class Squared_Exponential_Times_Sinusoidal_InflowField
 * 
 * \brief Implements an inflow field with the following expression: 
 *        \f$ Inflow(x,y) = M \left\{ A + \left( e^{S (r - r_{0})} \sin [S \pi (r - r_{0})] \right)^2 \right\} \f$ \n
 *
 * where r is the location of interest, \f$r_{0}\f$ is the reference point,
 * M is the magnitude of the function, A is a shift coefficient and S controls the
 * steepness and frequency. \n
 */
class Squared_Exponential_Times_Sinusoidal_InflowField: public InflowFieldBasicType{
public:
  
  //! Basic Constructor
  Squared_Exponential_Times_Sinusoidal_InflowField(void): ReferencePoint(0.0),
							  A(0.0), Magnitude(1.0),
							  SteepnessFrequency(1.0), Power(2)
  { 
    // Name the inflow field
    InflowFieldName = "Squared Exponential Sinusoidal, d=r-r0  [0:1], Inflow(x,y) = A + M*[exp(S*d)(sin(S*PI*d))^P], otherwise 0"; }
  
  //! Return value of the inflow field
  double EvaluateSolutionAt(const double &x, const double &y) const;

  //! Parse the input control parameters
  void Parse_Next_Input_Control_Parameter(AdvectDiffuse2D_Input_Parameters & IP, int & i_command);

  //! Print relevant parameters
  void Print_Info(std::ostream & out_file);
  
  //! Broadcast relevant parameters
  void Broadcast(void);
  
private:
  Vector2D ReferencePoint;	//!< 2D-space location used as reference for the field
  double A;			//!< Constant parameter of the inflow field
  double Magnitude, SteepnessFrequency; 	//!< Function magnitude and steepness-frequency control
  int Power;			//!< Function power
};

//! Return value of the inflow field
inline double Squared_Exponential_Times_Sinusoidal_InflowField::EvaluateSolutionAt(const double &x, const double &y) const {

  double R;
  double Func(1.0), TempFunc;
  int i;
  Vector2D PointOfInterest(x,y);

  // Calculate the distance R relative to the reference point
  R = PointOfInterest.abs() - ReferencePoint.abs();

  if (R >= 0 && R <= 1){
    // Calculate the sinusoidal
    TempFunc = sin(SteepnessFrequency * PI * R);

    // Calculate the power
    for (i = 1; i<= Power; ++i){
      Func *= TempFunc;
    }

    // Calculate main expression
    Func  *= exp(SteepnessFrequency*R);

    // Calculate final expression
    return A + Magnitude * Func;
  } else {
    return A;
  }
}

#endif
