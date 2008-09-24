/*!\file TestFunction_1D.cc
  \brief Source file to implement the test functions prototyped in TestFunctions_1D.h. */

/* Include required C++ libraries. */
#include <iostream>

/* Using std namespace functions */
using std::cout;

/* Include CFFC header files */
#include "TestFunctions_1D.h"
#include "../../Math.h"

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
    f =1.05*100;
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

  double f(0.0);
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

  double f(0.0);
  cout << "Theoretic integration Test_Example6 not defined!\n Please use Numeric Integration!!\n";
  return f;
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

  double f(0.0);
  cout << "Theoretic integration Test_Example8 is not defined!\n Please use Numeric Integration!!\n";
  return f;
}

double Test_Example9 (double x) {

  return sin(x*log(x));
}

double Test_Example10 (double x){
  return exp(x);
}

double Test_Example10_Integral (double x1, double x2){
  return exp(x2) - exp(x1);
}

double Test_Example11 (double x){
  return sqrt(x);
}

double Test_Example11_Integral (double x1, double x2){
  return (pow(x2,1.5) - pow(x1,1.5))/1.5;
}

double Test_Example12 (double x){
  return 1.0/(pow(x,4) + x*x + 0.9);
}

double Test_Example12_Integral (double x1, double x2){
return (.3096328464*atan(-.5565231379+1.174974040*x2)+.3096328464*atan(.5565231379+1.174974040*x2)-
	.3096328464*atan(-.5565231379+1.174974040*x1)-.3096328464*atan(.5565231379+1.174974040*x1)-
	.2781850613*log(9.486832981-9.472943556*x2+10.*pow(x2,2))+.2781850613*log(9.486832981+9.472943556*x2+10.*pow(x2,2))+
	.2781850613*log(9.486832981-9.472943556*x1+10.*pow(x1,2))-.2781850613*log(9.486832981+9.472943556*x1+10.*pow(x1,2)));
}

double Test_Example13 (double x){
  return (23.0/25.0)*cosh(x) - cos(x);
}

double Test_Example13_Integral (double x1, double x2){
  return (23.0/25.0)*(sinh(x2)-sinh(x1))-sin(x2)+sin(x1);
}

double Test_Example14 (double x){
  return sqrt(x*x*x);
}

double Test_Example14_Integral (double x1, double x2){
  return 0.4*(pow(x2,2.5) - pow(x1,2.5));
}

double Test_Example15 (double x){
  return 1.0/(1.0 + x);
}

double Test_Example15_Integral (double x1, double x2){
  return log(1.0+x2)-log(1.0+x1);
}

double Test_Example16 (double x){
  return 1.0/(1.0 + exp(x));
}

double Test_Example16_Integral (double x1, double x2){
  return log(exp(x2))-log(1.0+exp(x2))-log(exp(x1))+log(1.0+exp(x1));
}

double Test_Example17 (double x){
  return 1.0/(1.0 + pow(x,4));
}

double Test_Example17_Integral (double x1, double x2){
  double Sqrt2 = sqrt(2.0);
  return (0.25*Sqrt2)*(0.5*log(((x2*x2+Sqrt2*x2+1.0)*(-x1*x1+Sqrt2*x1-1.0))/((-x2*x2+Sqrt2*x2-1.0)*(x1*x1+Sqrt2*x1+1.0))) +
		       atan(Sqrt2*x2+1.0) - atan(Sqrt2*x1+1.0) + atan(Sqrt2*x2-1.0) - atan(Sqrt2*x1-1.0));
}

double Test_Example18 (double x){
  return 1.234*pow(x,18) + 0.5*pow(x,15) - 0.25*pow(x,8) + 0.0000234*x*x - 0.000001*x + 0.001;
}

double Test_Example18_Integral (double x1, double x2){
  double f(0.0);
  cout << "Theoretic integration Test_Example8 is not defined!\n Please use Numeric Integration!!\n";
  return f;
}

double Test_Example19 (const double & x){
  return exp(-2.0*(x+5.0)) - log(x+1.0)*sin(x);
}
