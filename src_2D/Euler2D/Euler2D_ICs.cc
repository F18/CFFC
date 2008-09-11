/*!\file Euler2D_ICs.cc
   \brief Source file initializing/implementing member variables/functions define in Euler2D_ICs.h file. */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "Euler2D_ICs.h"

double SinVariationInXDir(const double x, const double y){
  if (x<-100 || x>100){
    return 2.0;
  } else {
    return 2.0 + std::sin((ConvertDomainToMinusOneOne(-100,100,x)+1)*PI);
  }
}


double SinVariationInYDir(const double x, const double y){
  if (y<-100 || y>100){
    return 2.0;
  } else {
    return 2.0 + std::sin((ConvertDomainToMinusOneOne(-100,100,y)+1)*PI);
  }
}

double MultipleSinVariationInXDir(const double x, const double y){

  return (1.0 + 0.5*std::cos(PI*ConvertDomainToMinusOneOne(-100,100,x))*
	  std::sin(5*ConvertDomainToMinusOneOne(-100,100,x)*PI) );

}

double MultipleSinVariationInYDir(const double x, const double y){

  return (1.0 + 0.5*std::cos(PI*ConvertDomainToMinusOneOne(-100,100,y))*
	  std::sin(5*ConvertDomainToMinusOneOne(-100,100,y)*PI) );
}


double Complex_2D_Waves(const double x, const double y){

  // Gaussian function
  double A(100), A_1(350);
  double x0(-50.0), y0(-50.0), x0_1(40.0), y0_1(30.0);
  double theta(PI/5), theta_1(PI/2);
  double sigma_x(5.7);
  double sigma_y(10.0), sigma_1(13.0);
  double a,b,c, tmp1, tmp2, Z, Z_1;

  // compute first Gaussian function
  tmp1 = (cos(theta)/sigma_x);
  tmp2 = (sin(theta)/sigma_y);
  a = tmp1 * tmp1 + tmp2* tmp2;
  
  b = -sin(2*theta)/(sigma_x*sigma_x) + sin(2*theta)/(sigma_y*sigma_y);

  tmp1 = (sin(theta)/sigma_x);
  tmp2 = (cos(theta)/sigma_y);
  c = tmp1 * tmp1 + tmp2 * tmp2;

  Z = A*exp( - (a*(x-x0)*(x-x0) + b*(x-x0)*(y-y0) + c*(y-y0)*(y-y0)));


  // compute second Gaussian function
  tmp1 = (cos(theta_1)/sigma_1);
  tmp2 = (sin(theta_1)/sigma_1);
  a = tmp1 * tmp1 + tmp2* tmp2;
  
  tmp1 = (sin(theta_1)/sigma_1);
  tmp2 = (cos(theta_1)/sigma_1);
  c = tmp1 * tmp1 + tmp2 * tmp2;

  Z_1 = A_1*exp( - (a*(x-x0_1)*(x-x0_1) + c*(y-y0_1)*(y-y0_1)));


  return (50.0 + Z + Z_1 );

  //   return (50.0 + Z + Z_1 + 15.0*std::cos(3*ConvertDomainToMinusOneOne(-100,100,x))*
  // 	  std::sin(3*ConvertDomainToMinusOneOne(-50,50,y)*PI)*
  // 	  std::sin(5*ConvertDomainToMinusOneOne(-100,100,x)*PI) );
}

double Helper_Abgrall_2D_Function(const double r)
{
  if (r <= -0.1e1 / 0.3e1){
    return (-r * std::sin(0.3e1 / 0.2e1 * PI * r * r));
  } else if (std::fabs(r) < 0.1e1 / 0.3e1){
    return (std::fabs(std::sin(0.2e1 * PI * r)));
  } else if (0.1e1 / 0.3e1 <= r){
    return (0.2e1 * r - 0.1e1 + std::sin(0.3e1 * PI * r) / 0.6e1);
  } else {
    return 0.0e0;
  }
}

double Abgrall_2D_Function(const double x, const double y)
{
  if (x <= std::cos(PI * y) / 0.2e1){
    return (2.0 + Helper_Abgrall_2D_Function(x - 1.0/std::tan(std::sqrt(PI / 0.2e1)) * y));
  } else if (std::cos(PI * y) / 0.2e1 < x) {
    return (2.0 + Helper_Abgrall_2D_Function(x + 1.0/std::tan(std::sqrt(PI / 0.2e1)) * y) +
	    std::cos(0.2e1 * PI * y));
  } else {
    return 0.0e0;
  }
}

double SinExponentialVariationInXDir(const double x, const double y){
  return 500.0 + 1.0e-5*(exp(-4*x)*sin(5*x));
}

double SinExponentialVariationInYDir(const double x, const double y){
  return 500.0 + 1.0e-5*(exp(-4*y)*sin(5*y));
}

double SinExponentialVariationRotated(const double x, const double y){
  std::cout << "Not implementated yet\n";
  return 0.0;
}


double CosineHill(const double x, const double y){

  static Vector2D Centroid(-20,-20);
  static double ConeRadius(15);

  Vector2D CurrentNode(x,y);

  double Distance;

  Distance = abs(CurrentNode - Centroid);

  if (Distance <= ConeRadius){	// point inside or on the edge of the base circle
    return 2.0 + 0.5*(1.0 + cos(PI/ConeRadius*Distance));
  }

  return 2.0; 			// point outside the base circle
}

double Translated_CosineHill(const double x, const double y){

  std::cout << "Centroid of the translated function not defined yet!\n";

  static Vector2D Centroid(-20,-20);
  static double ConeRadius(15);

  Vector2D CurrentNode(x,y);

  double Distance;

  Distance = abs(CurrentNode - Centroid);

  if (Distance <= ConeRadius){	// point inside or on the edge of the base circle
    return 2.0 + 0.5*(1.0 + cos(PI/ConeRadius*Distance));
  }

  return 2.0; 			// point outside the base circle
}

double Polynomial_Function(const double x, const double y){

  return 200.0 + (x-0.25)*(x-0.5)*(x+0.25)*(y-0.24)*(y-1.5)*(y+0.25);

}

// Initialization of static member variables in Translated_Solutions class
Vector2D Translated_Solutions::InitialCentroid = Vector2D(0.0);
Vector2D Translated_Solutions::TranslatedCentroid = Vector2D(0.0);
double Translated_Solutions::Magnitude = 1.0;
double Translated_Solutions::Steepness = 1.0;
double Translated_Solutions::ConeRadius = 1.0;
double Translated_Solutions::Distance = 0.0;

