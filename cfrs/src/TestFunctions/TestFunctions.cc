/* FunctionTest.cc: Source file defining
                    the forms of the test functions for reconstruction 1D and 2D*/

/* Include header files. */

#include "TestFunctions.h"

/*******************************************************
 *   !!!!!!!!!!!!!   ATTENTION !!!!!!!!!!!!!           *
 *                                                     *
 *   x1 -- the lower limit                             *
 *   x2 -- the upper limit                             *
 ******************************************************/


/********************************************************
 * Function Test_Default1D :                            *
 *          f(x) = x^6 - 4x^4 + 2x                      *
 *******************************************************/

double Test_Default1D (double x){

  double f;
  f = (pow(x,6)-4*pow(x,4)+2*x);
  return f;
}
double Test_Default1D_Integral (double x1, double x2) {

  double f;
  f = (1.0/7.0*(pow(x2,7) - pow(x1,7)) - 4.0/5.0*(pow(x2,5)-pow(x1,5)) + 
      (x2*x2 - x1*x1));
  return f;
}
double Test_Default1D_FirstDeriv (double x){

  double f;
  f = 6*pow(x,5)-16*pow(x,3)+2;
  return f;
}
double Test_Default1D_SecondDeriv (double x){

  double f;
  f = 30*pow(x,4)-48*pow(x,2);
  return f;
}
double Test_Default1D_ThirdDeriv (double x){

  double f;
  f = 120*pow(x,3)-96*x;
  return f;
}
double Test_Default1D_FourthDeriv (double x){

  double f;
  f = 360*pow(x,2)-96;
  return f;
}


/********************************************************
 * Function Test_Example1:                              *
 *          f(x) = sin(2x)+2*cos(x)                     *
 *******************************************************/

double Test_Example1 (double x) {

  double f;
  f = (sin(2*x)+2*cos(x));
  return f;
}
double Test_Example1_Integral (double x1, double x2) {

  double f;
  f = -0.5*(cos(2*x2)-cos(2*x1))+2*(sin(x2)-sin(x1));
  return f;
}
double Test_Example1_FirstDeriv (double x) {

  double f;
  f = 2*cos(2*x)-2*sin(x);
  return f;
}
double Test_Example1_SecondDeriv (double x) {

  double f;
  f = -4*sin(2*x)-2*cos(x);
  return f;
}
double Test_Example1_ThirdDeriv (double x) {

  double f;
  f = -8*cos(2*x)+2*sin(x);
  return f;
}
double Test_Example1_FourthDeriv (double x) {

  double f;
  f = 16*sin(2*x)+2*cos(x);
  return f;
}

/********************************************************
 * Function Test_Example2:                              *
 *          f(x) = exp(-4*x)*sin(5*x)                   *
 *******************************************************/

double Test_Example2 (double x) {

  double f;
  f = exp(-4*x)*sin(5*x);
  return f;
}

double Test_Example2_g (double x) {
  // defined function for computing the integral of Test_Example2
  double g;
  g = exp(-4*x)*cos(5*x);
  return g;
}

double Test_Example2_Integral (double x1, double x2) {
  double f;
  f = -5.0/41.0*(Test_Example2_g(x2)-Test_Example2_g(x1)) - 
    4.0/41.0*(Test_Example2(x2)-Test_Example2(x1));
  return f;
}

double Test_Example2_FirstDeriv (double x) {
  double f;
  f = -20*Test_Example2_g(x);
  return f;
}

double Test_Example2_SecondDeriv (double x) {
  double f;
  f = -400*Test_Example2(x);
  return f;
}

double Test_Example2_ThirdDeriv (double x) {
  double f;
  f = 8000*Test_Example2_g(x);
  return f;
}

double Test_Example2_FourthDeriv (double x) {
  double f;
  f = 160000*Test_Example2(x);
  return f;
}

/********************************************************
 * Function Test_Example3:                              *
 *          f(x) = step function at x = 2               *
 *******************************************************/

static double StepPoint = 2.0;
static double TolStepPoint = 1.0e-10;

double Test_Example3 (double x) {
  double f;
  if(x < (StepPoint-TolStepPoint))
    f =1.05e8;
  else 
    f = 1.0;
  return f;
}

double Test_Example3_Integral (double x1, double x2) {
  double f;
  if (x2 <= (StepPoint-TolStepPoint))
    f = (x2-x1)*Test_Example3(x2);
  else if (x1<=(StepPoint-TolStepPoint))
    f = (StepPoint-x1)*Test_Example3(x1) + (x2-StepPoint)*Test_Example3(x2);
  else
    f = (x2 - x1)*Test_Example3(x1);
  return f;
}

/********************************************************
 * Function Test_Example4:                              *
 *          f(x) = 2.0 + sin((x+1)* PI)                 *
 *******************************************************/

double Test_Example4 (double x) {

  double f;
  f = 2.0 + sin((x+1)* PI);
  return f;
}

double Test_Example4_Integrand (double x) {

  double f;
  f = 2*x + cos(PI*x)/PI;

  return f;
}

double Test_Example4_Integral (double x1, double x2) {

  double f;
  f = Test_Example4_Integrand(x2) - Test_Example4_Integrand(x1);

  return f;
}

/********************************************************
 * Function Test_Example5:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example5 (double x) {

  double f;
  double val = 1.0/3.0;  

  if (x <= -val)
    f = -x*sin(3.0*PI*x*x/2.0);
  else if (x < val)
    f = fabs(sin(2*PI*x));
  else
    f = 2*x-1+1.0/6.0*sin(3.0*PI*x);
 
  return f;
}

double Test_Example5_Integral (double x1, double x2) {

  double f;
  f = x1; f = x2; f=0.0;		// that's for the compiler!
  cout << "Theoretic integration Test_Example5 not defined!\n Please use Numeric Integration!!\n";
  return f;
}

/********************************************************
 * Function Test_Example6:                              *
 *          f(x) = 0.345*x + 0.5                        *
 *******************************************************/

double Test_Example6 (double x) {

  double f;
  f = 0.345*x + 0.5;
  return f;
}

double Test_Example6_Integral (double x1, double x2) {

  double f;
  cout << "Theoretic integration Test_Example6 not defined!\n Please use Numeric Integration!!\n";
  return 0;
}

/********************************************************
 * Function Test_Example7:                              *
 *          f(x) = sin(PI*x)                            *
 *******************************************************/

double Test_Example7 (double x) {

  double f;
  f = sin(PI*x);
  return f;
}

double Test_Example7_Integral (double x1, double x2) {

  double f;
  f = -1.0/PI * (cos(PI*x2) - cos(PI*x1));
  return f;
}

/********************************************************
 * Function Test_Example8:                              *
 *          f(x) = 2 steps of different intensities     *
 *******************************************************/

double Test_Example8 (double x) {

  double f;

  if(x < 1){
    f =1.0;
  } else if (x < 2) {
    f = 1.4*7.8e3;
  }
  else {
    f = 1.3e5;
  }
  return f;
}

double Test_Example8_Integral (double x1, double x2) {

  double f;
  f = x1; f = x2; f=0.0;		// that's for the compiler!
  cout << "Theoretic integration Test_Example8 is not defined!\n Please use Numeric Integration!!\n";
  return f;
}

/********************************************************
 * Function Test_Example9:                              *
 *          f(x) =                    *
 *******************************************************/

double Test_Example9 (double x){

  double f;
  double step = 1.05e4;
  double xCoord;

  if(x < StepPoint - 0.1){
    xCoord = x - (StepPoint - 0.1) + 0.5*PI;
    f = step + 100*(sin(2*xCoord)+2*cos(xCoord));
  } else if(x < (StepPoint-TolStepPoint))
    f =step;
  else 
    f = 1.0;

  return f;
}

double Test_Example9_Integral (double x1, double x2){
  double f;
  f = x1; f = x2; f=0.0;		// that's for the compiler!
  cout << "Theoretic integration Test_Example9 is not defined!\n Please use Numeric Integration!!\n";
  return f;
}


/********************************************************
 * Function Test_Example10: Jiang function              *
 *          f(x) =                                      *
 *******************************************************/
double Test_Example10 (double x){
  /*Obs: The function is shifted on "X" and "Y"
    X-shift -> in order to position the singular points inside the interval
    Y-shift -> in order to use this data distribution as initial condition for density in Euler1D
  */

  return 6.0 + JIANG_Function(x-0.5);
}

double Test_Example10_Integral (double x1, double x2){
  double f;
  f = x1; f = x2; f=0.0;		// that's for the compiler!
  cout << "Theoretic integration Test_Example10 is not defined!\n Please use Numeric Integration!!\n";
  return f;
}

double JIANG_Function (double x){

  double point1, point2, point3, point4, point5;
  point1 = -1.0;
  point2 = - 1.0/3.0;
  point3 =  0.0;
  point4 = 1.0/3.0;
  point5 = 1.0;

  if ( x >= -1.5 && x < point1){
    x += 2.0; // periodic function
  }

  double solution = -(sqrt(3.0)/2.0 + 4.5 + 2.0*PI/3.0)*(x+1);

  if( point1 <= x && x < point2 ) {
    solution += 2.0*cos(3.0*PI*x*x/2.0) - sqrt(3.0);
    
  } else if( point2 <= x && x < point3 ) {
    solution += 3.0/2.0 + 3.0*cos(2.0*PI*x);

  } else if( point3 <= x && x < point4 ){
    solution += 15.0/2.0 - 3.0*cos(2.0*PI*x);

  } else if( point4 <= x && x <= point5 ) {
    solution += (28.0 + 4*PI + cos(3.0*PI*x))/3.0 + 6*PI*x*(x-1);

  } else {
    solution = 0.0;
  } //endif

  return solution;
}
