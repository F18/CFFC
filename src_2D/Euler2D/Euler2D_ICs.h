/*!\file Euler2D_ICs.h
   \brief Header file defining 2D Euler initial conditions defined analytically. */

#ifndef _EULER2D_ICS_INCLUDED
#define _EULER2D_ICS_INCLUDED

/* Include required C++ libraries. */
#include <cmath>
#include <iostream>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"
#include "../Math/Vector2D.h"

double SinVariationInXDir(const double x, const double y);

double SinVariationInYDir(const double x, const double y);

double MultipleSinVariationInXDir(const double x, const double y);

double MultipleSinVariationInYDir(const double x, const double y);

double Complex_2D_Waves(const double x, const double y);

double Abgrall_2D_Function(const double x, const double y);

double Helper_Abgrall_2D_Function(const double r);

double SinExponentialVariationInXDir(const double x, const double y);

double SinExponentialVariationInYDir(const double x, const double y);

double SinExponentialVariationRotated(const double x, const double y);

double CosineHill(const double x, const double y);

double Translated_CosineHill(const double x, const double y);

double Polynomial_Function(const double x, const double y);

class Translated_Solutions{
  
public:

  // Points defining the translation of the function
  static Vector2D InitialCentroid;
  static Vector2D TranslatedCentroid;

  // Control parameters
  static double Magnitude, Steepness;

  // Other parameters
  static double ConeRadius;
  static double Distance;

  // ============ INITIAL CONDITIONS ===============

  // Translation of function "1 - tanh(DistanceToCentroid)^2"
  static double HyperTangentIC(double x, double y);
  static double Translated_HyperTangent(double x, double y);

 
  
private:
  Translated_Solutions(){ };
  
  // Helper functions
  static double HyperTangent(const double & DistanceToCentroid);

};

// Definitions of member functions of Translated_Solutions
inline double Translated_Solutions::HyperTangentIC(double x, double y){

  Distance = abs(InitialCentroid - Vector2D(x,y));

  return HyperTangent(Distance);
}

inline double Translated_Solutions::Translated_HyperTangent(double x, double y){

  Distance = abs(TranslatedCentroid - Vector2D(x,y));

  return HyperTangent(Distance);
}

inline double Translated_Solutions::HyperTangent(const double & DistanceToCentroid){

  // Function: "1 - tanh(x)^2" modulated for magnitude and steepness

  double Func;

  Func = tanh(Steepness * DistanceToCentroid);
  Func *= Func; 		// square hyperbolic tangent
  Func = 1.0 - Func;

  return 1.0 + Magnitude*Func;
}


#endif
