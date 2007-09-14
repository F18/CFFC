/* FunctionTest.cc: Source file defining
                    the forms of the test functions for reconstruction 2D*/

/* Include header files. */

#ifndef _TESTFUNCTIONS_2D_INCLUDED
#include "TestFunctions_2D.h"
#endif // _TESTFUNCTIONS_2D_INCLUDED

#include <cmath>
#include <cassert>

/********************************************************
 * Function Test_Default1D :                            *
 *          f(x) = x^6 - 4x^4 + 2x                      *
 *******************************************************/

double Test_Default2D (double x, double y){

  double f;
  f = x;			// that's just for the compiler!!!
  //  f = pow(x,6)-4*pow(x,4)+2*x;
  f = pow(y,6)-4*pow(y,4)+2*y;
  return f;
}

//this expression's only for cartesian mesh !!!!!!!!!!!!!!!!!
double Test_Default2D_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double f;
  f = 1.0/7.0*(pow(y2,7) - pow(y1,7)) - 4.0/5.0*(pow(y2,5)-pow(y1,5)) + 
    (y2*y2 - y1*y1);
  f *= (x2 - x1);
  return f;
}

/********************************************************
 * Function Test_Example1:                              *
 *          f(x,y) = x*x + y*y + x*y                    *
 *******************************************************/

double Test_Example1 (double x, double y) {
  double f;

  f = x*x + y*y +x*y;
  return f;
}

double Test_Example1_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double f;
  //f = 1.0/3.0*((pow(x2,3) - pow(x1,3))*(y2-y1) +(x2-x1)*(pow(y2,3) - pow(y1,3)));
  // for x^2+y^2 

  f = 1.0/3.0*(x2-x1)*(pow(y2,3)-pow(y1,3)) + 0.25*(x2*x2-x1*x1)*(y2*y2-y1*y1) + 
    1.0/3.0*(pow(x2,3)- pow(x1,3))*(y2-y1);

  return f;
}

/**********************************************************************
 * Function Test_Example2: 6th order polynomial function              *
 *          f(x) = (x-0.25)*(x-0.5)*(x+0.25)*(y-0.24)*(y-1.5)*(y+0.25)*
 *********************************************************************/

double Test_Example2 (double x, double y) {

  double f;

  f = (x-0.25)*(x-0.5)*(x+0.25)*(y-0.24)*(y-1.5)*(y+0.25);
  return f;
}

double Test_Example2_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double func;
  double a, b, c, d, e, f, g,h, j, k, l, m, n;

  a = 0.25;
  b = 0.1666666666;
  c = 0.03125;
  d = 0.3333333333;
  e = 0.3725;
  f = 0.2483333333;
  g = 0.0465625;
  h = 0.01875;
  j = 0.0125;
  k = 0.00234375;
  l = 0.0225;
  m = 0.015;
  n = 0.0028125;

  func = a*(a*pow(x2,4)-a*pow(x1,4)-b*pow(x2,3)+b*pow(x1,3)-c*pow(x2,2)+
	    c*pow(x1,2)+c*x2-c*x1)*(pow(y2,4)-pow(y1,4))+
         d*(-e*pow(x2,4)+e*pow(x1,4)+f*pow(x2,3)-f*pow(x1,3)+g*pow(x2,2)-
	    g*pow(x1,2)-g*x2+g*x1)*(pow(y2,3)-pow(y1,3))+
         0.5*(-h*pow(x2,4)+h*pow(x1,4)+j*pow(x2,3)-j*pow(x1,3)+k*pow(x2,2)-
	    k*pow(x1,2)-k*x2+k*x1)*(pow(y2,2)-pow(y1,2))+
        (l*pow(x2,4)-l*pow(x1,4)-m*pow(x2,3)+m*pow(x1,3)-n*pow(x2,2)+
	 n*pow(x1,2)+n*x2-n*x1)*(y2-y1);

  return func;
}

/*********************************************************
 * Function Test_Example3: 5th order polynomial function *
 *            f(x) = x^4*y+y^4+8*x^2*y-6*y^3             *
// *          f(x) = (x-0.5)*(y-0.25)+(y+2)^3            *
 ********************************************************/

double Test_Example3 (double x, double y) {

  double f;
  //  f = (x-0.5)*(y-0.25)+pow(y+2,3);
  f = pow(x,4)*y+ pow(y,4)+8*x*x*y-6*pow(y,3);
  return f;
}

double Test_Example3_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double func;

//    func = 0.25*((x2-x1)*(pow(y2,4)-pow(y1,4)) + (x2*x2-x1*x1)*(y2*y2-y1*y1))+
//      2.0*(x2-x1)*(pow(y2,3)-pow(y1,3)) + 5.75*(x2-x1)*(y2*y2-y1*y1)+
//      0.125*(x2*x2-x1*x1)*(-y2+y1)+ 8.125*(x2-x1)*(y2-y1);

  func = 0.2*(x2-x1)*(pow(y2,5)-pow(y1,5)) + 1.5*(-x2+x1)*(pow(y2,4)-pow(y1,4))+
    0.5*(0.2*pow(x2,5)-0.2*pow(x1,5) + 8.0/3.0*pow(x2,3) - 8.0/3.0*pow(x1,3))*
    (y2*y2 - y1*y1);

  return func;
}

/********************************************************
 * Function Test_Example4:                              *
 *          f(x) = x*(x-0.25)*(x+0.25)*pow(y,3)         *
 *******************************************************/

double Test_Example4 (double x, double y) {

  double f;
  f = x*(x-0.25)*(x+0.25)*pow(y,3);
  return f;
}

double Test_Example4_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double f;
  f = 0.0078125*(-8*pow(x2,4)+8*pow(x1,4)+x2*x2-x1*x1)*(-pow(y2,4)+pow(y1,4));
  return f;
}

/********************************************************
 * Function Test_Example5:                              *
 *          f(x) = step function at x=0.5               *
 *          (only in one direction)                     *
 *******************************************************/

double Test_Example5 (double x, double y) {
  double f;
  double a = 0.5;
  f = x; 			// that's just for the compiler
  if(y <= a)
    f =1.0e5;
  else 
    f = 1.0;
  return f;
}

double Test_Example5_Integral (double x1, double x2, double y1, double y2) {
  assert((x1<x2)&&(y1<y2));
  double tol = 0.1e-6;
  double f;

  double a = 0.5;

  if (y2 <= (a-tol))
    f = (y2-y1)*Test_Example5(x2,y2);
  else if ((y1<=a)&&(y2>a))
    f = (a-y1)*Test_Example5(x1,y1) + (y2-a)*Test_Example5(x2,y2);
  else
    f = (y2 - y1)*Test_Example5(x1,y1);
  f *= (x2 - x1);
  return f;

}

/********************************************************
 * Function Test_Example6:                              *
 *          f(x) = step function in 2D                  *
 *******************************************************/

double Test_Example6 (double x, double y) {

  double f;
  double x_step = 0.5;
  double y_step = 0.3;

  if (x<x_step){
    f = 1.023;
  } else if (y<y_step)
    f = 1.023;
  else
    f = 10.423*1.0e3;
  return f;
}

double Test_Example6_Integral (double x1, double x2, double y1, double y2) {
  double f;
  double x_step = 0.5;
  double y_step = 0.1;

  assert ((x1<x2)&&(y1<y2));
  if ((x2<x_step)||(y2<y_step)||((x1>=x_step)&&(y1>=y_step))){
    f = (x2-x1)*(y2-y1)*Test_Example6(x1,y1);
  } else if ((x1<x_step)&&(x2>=x_step)){
    if (y1 < y_step){
      f = (x_step-x1)*(y2-y1)*Test_Example6(x1,y1);
      f += (x2-x_step)*(y_step-y1)*Test_Example6(x2,y1);
      f += (x2-x_step)*(y2-y_step)*Test_Example6(x2,y2); 
    }
    else {
      f = (x_step-x1)*(y2-y1)*Test_Example6(x1,y1);
      f += (x2-x_step)*(y2-y1)*Test_Example6(x2,y2);
    }
  } else{
    f = (y2-y_step)*Test_Example6(x2,y2)+(y_step-y1)*Test_Example6(x1,y1);
    f *= (x2-x1);
  }
  return f;
}

/********************************************************
 * Function Test_Example7:                              *
 *          f(x) = cos(PI*x^2+4*PI*y)                   *
 *******************************************************/

double Test_Example7 (double x, double y) {

  double f;
  f = cos(PI*x*x+4*PI*y);
  return f;
}

double Test_Example7_Integral (double x1, double x2, double y1, double y2) {

  assert ((x1<x2)&&(y1<y2));
  double f;
  double FresnelC_x2, FresnelC_x1, FresnelS_x2, FresnelS_x1;
  double Sqrt2 = sqrt(double(2));
  double val;

  frenel(Sqrt2*x2, &FresnelS_x2, &FresnelC_x2);
  frenel(Sqrt2*x1, &FresnelS_x1, &FresnelC_x1);

  val = 2*PI*(y2+y1);
  f = cos(val)*(FresnelC_x2 - FresnelC_x1)-sin(val)*(FresnelS_x2-FresnelS_x1);
  f *= Sqrt2*sin(2*PI*(y2-y1))/(4*PI);
  return f;
}

/********************************************************
 * Function Test_Example8 (Abgral Function):            *
 *          f(x) =                    *
 *******************************************************/

double Test_Example8_f (double r)
{
  double cgret;
  if (r <= -0.1e1 / 0.3e1)
    cgret = -r * sin(0.3e1 / 0.2e1 * PI * r * r);
  else if (fabs(r) < 0.1e1 / 0.3e1)
    cgret = fabs(sin(0.2e1 * PI * r));
  else if (0.1e1 / 0.3e1 <= r)
    cgret = 0.2e1 * r - 0.1e1 + sin(0.3e1 * PI * r) / 0.6e1;
  else
    cgret = 0.0e0;
  return(cgret);
}

double Test_Example8 (double x, double y)
{
  double cgret;
  if (x <= cos(PI * y) / 0.2e1)
    cgret = Test_Example8_f(x - 1.0/tan(sqrt(PI / 0.2e1)) * y);
  else if (cos(PI * y) / 0.2e1 < x)
    cgret = Test_Example8_f(x + 1.0/tan(sqrt(PI / 0.2e1)) * y) +
      cos(0.2e1 * PI * y);
  else
    cgret = 0.0e0;
  return(cgret);
}

double Test_Example8_Integral (double x1, double x2, double y1, double y2) {

  double f;

  // that's just for the compiler
  f = x1;
  f = x2;
  f = y1;
  f = y2;
  f = 0.0;

  return f;
}

double Test_Example9 (double x, double y){

  double f;
  double a = -0.75/0.7;
  double b = 0.75;
  double c = 1.2/3.5;
  double d = -0.86;
  double low_value = 2.3;
  double high_value = 1000.0;
  
  if ((y <= a*x+b)||(y <= c*x+d))
    f = low_value;
  else
    f = high_value;

  return f;
}

/********************************************************
 * Function Test_Example10:                             *
 *          f(x,y) = exp((x+y)/(x-y))                   *
 *******************************************************/

double Test_Example10(double x, double y){

  return exp((x+y)/(x-y));
}

/********************************************************
 * Function Test_Example11:                             *
 *          f(x,y) = 1/50*exp(-x/10)*exp(-y/5)          *
 *******************************************************/

double Test_Example11(double x, double y){

  return  1.0/50.0*exp(-x/10)*exp(-y/5);
}


/********************************************************
 * Function Test_Example12:                             *
 *          f(x,y) = step function                      *
 *******************************************************/

double Test_Example12(double x, double y){
  double f;
  // that's just for the compiler
  f = y;

  double a = 0.0;
  if(x > a)
    f =1.0e5;
  else 
    f = 1.0;
  return f;
}

/********************************************************
 * Function Test_Example12:                             *
 *          f(x,y) = x^3*y^5                            *
 *******************************************************/

double Test_Example13(double x, double y){

  return pow(x,3)*pow(y,5);
}

