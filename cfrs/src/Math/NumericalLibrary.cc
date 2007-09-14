/* NumericalLibrary.cc provides useful implementation for different numerical methods*/

#include "NumericalLibrary.h"

void frenel(double x, double *s, double *c)
{
  /* returns the sin (*s) and cos (*c) values for the Frenel function applied to "x"*/
  double EPS=6.0e-8;
  int MAXIT=100;
  double FPMIN=1.0e-30;
  double XMIN=1.5;
  double PIBY2 (PI/2.0);
  int TRUE=1;

  int k,n,odd;
  float a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
  complex<double> b, cc, d, h, del, cs, a_compl;
  
  ax=fabs(x);
  if (ax < sqrt(FPMIN)) {

    *s=0.0;
    *c=ax;
  } else if (ax <= XMIN) {
    sum=sums=0.0;
    sumc=ax;
    sign=1.0;
    fact=PIBY2*ax*ax;
    odd=TRUE;
    term=ax;
    n=3;
    for (k=1;k<=MAXIT;k++) {
      term *= fact/k;
      sum += sign*term/n;
      test=fabs(sum)*EPS;
      if (odd) {
	sign = -sign;
	sums=sum;
	sum=sumc;
      } else {
	sumc=sum;
	sum=sums;
      }
      if (term < test) break;
      odd=!odd;
      n += 2;
    }
    if (k > MAXIT) cout << "series failed in frenel" << endl;
    *s=sums;
    *c=sumc;
  } else {
    pix2 = PI*ax*ax;
    b = complex<double>(1.0,-pix2);
    cc = 1.0/FPMIN;
    d = h = 1.0/b;
    n=-1;
    for (k=2; k<=MAXIT; k+1){
      n +=2;
      a = -n*(n+1);
      b += 4.0;
      a_compl = complex<double>(a,0.0);
      d = 1.0/(a_compl*d+b);
      cc = b + a_compl/cc;
      del = cc*d;
      h *= del;
      if (fabs(del.real()-1.0) + fabs(del.imag()) < EPS)
	break;
    }
    if (k > MAXIT) 
      cout << "cf failed in frenel" << endl;

    h *=complex<double>(ax,-ax);
    cs = complex<double>(0.5,0.5)*
      (complex<double>(1.0,0.0)-h*complex<double>(cos(0.5*pix2),sin(0.5*pix2)));
    *c = cs.real();
    *s = cs.imag();
  }
  if (x < 0.0) {
    *c = -(*c);
    *s = -(*s);
  }
  return;
}

double qgaus10(const FunctionType1D func, const double a, const double b){
     /* returns the integral of the function "func" between a and b,
	by ten-point Gauss-Legendre integration: the function is evaluated exactly
	ten times at interior points in the range of integration.
	The integral is exact for polynomials up to order of 19. 
	Implementation based on the subroutine from Numerical Recipes*/

	int j;
	double xr,xm,dx,s;
	static double x[]={0.0,0.148874338981631,0.433395394129247,
			   0.67940956829902444,0.86506336668898454,0.97390652851717174};
	static double w[]={0.0,0.29552422471475281,0.26926671930999635,
			   0.21908636251598207, 0.14945134915058053, 0.066671344308688082};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=1;j<=5;j++) {
		dx=xr*x[j];
		s += w[j]*(func(xm+dx)+func(xm-dx));
	}
	return s *= xr;
}

double qgaus5(const SuperFunctionType1D func, const FunctionType1D SubFunc1,
	      const FunctionType1D SubFunc2, const double a, const double b){
     /* returns the integral of the function "func" between a and b,
	by ten-point Gauss-Legendre integration: the function is evaluated exactly
	ten times at interior points in the range of integration.
	The integral is exact for polynomials up to order of 19. 
	Implementation based on the subroutine from Numerical Recipes*/

	int j;
	double xr,xm,dx,sum;
	static double x[]={0.0, 0.53846931010568311, 0.90617984593866396};
	static double w[]={0.56888888888888889, 0.47862867049936647, 0.23692688505618917};

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	sum=w[0]*(func(xm,SubFunc1,SubFunc2));
	for (j=1;j<=2;j++) {
		dx=xr*x[j];
		sum += w[j]*(func(xm+dx,SubFunc1,SubFunc2)+func(xm-dx,SubFunc1,SubFunc2));
	}
	sum *= xr;
	return sum;
}

double AdaptiveGaussianQuadrature(const SuperFunctionType1D func, const FunctionType1D SubFunc1,
				  const FunctionType1D SubFunc2, const double a,const double b,
				  int digits)
{
  double EPS = 0.05*pow(10.0,1.0-digits); // accuracy: --> based on precision 
                                         // i.e exact number of digits 
  double g5;                            // the value of the integral obtained
                                         // with qgaus5
  double RERR=1;	        	 // Relative Error
  long double G=0;			// Integral value
  int n;
  double q;
  long double LengthInterval = fabs(b - a); 
  long double StartPoint = a;
  long double EndPoint = b;
  long double SubDivisionPoint;
  double Division = 0.5;

  try
    {
      if (0.005*LengthInterval + 1 == 1)
	throw too_short_interval();
      g5 = qgaus5 (func,SubFunc1,SubFunc2,StartPoint,EndPoint);
      SubDivisionPoint = StartPoint + Division*(EndPoint - StartPoint);
      G = qgaus5(func,SubFunc1,SubFunc2,StartPoint,SubDivisionPoint)+ 
	       qgaus5(func,SubFunc1,SubFunc2,SubDivisionPoint,EndPoint);
      RERR = fabs(G - g5)/(1 + fabs(G));
       if (RERR < EPS){
	return G;
      }
      else {
	G = AdaptiveGaussianQuadrature(func,SubFunc1,SubFunc2,StartPoint,SubDivisionPoint,digits) + 
	  AdaptiveGaussianQuadrature(func,SubFunc1,SubFunc2,SubDivisionPoint,EndPoint,digits);
	return G;
      }
    }
  
  catch (const maximum_exceeded &)
    {
      std::string msg("Maximum iterations exceeded. The precision couldn't be obtained!");
      std::cerr << msg << std::endl;
      return 0;
    }

  catch (const machine_accuracy_exceeded &)
    {
      std::string msg("The machine accuracy exceeded! The precision cannot be obtained!");
      std::cerr << msg << std::endl;
      return 0;
    }

  catch (const too_short_interval &)
    {
      std::string msg("Too short interval! The precision cannot be obtained!");
      std::cerr << msg << std::endl;
      exit(0);
    }

  catch (...)
    {
      std::cerr << "Undefined error!" << std::endl;
    }

  return G;
}

static double xsav,yy1,yy2;
static FunctionType2D FuncPtr;

double quad2d(const FunctionType2D func, const double a, const double b, const double c, const double d){
/*returns the integral of a user-supplied func over a two-dimensional 
  rectangular region, specified by the limits x1, x2, y1, y2.
  Integration is performed by calling qgaus10 recursively.*/

  FuncPtr = func;
  yy1 = c; yy2 = d;
  return qgaus10(f1,a,b);
}
 
double f1(const double x){ 
  /*        /y2
    f1 = Int|  func(x,y)dy
            /y1        
  */

  xsav=x;
  return qgaus10(f2,yy1,yy2);
}

double f2(const double y){
  // The integrand f(x,y) evaluated at fixed x and y
  return FuncPtr(xsav,y);
}

static int NbDigitsAGQ;

double quad2dAdaptiveGaussianQuadrature(const FunctionType2D func, const double a,
					const double b, const double c, const double d, int digits){
/*returns the integral of a user-supplied func over a two-dimensional 
  rectangular region, specified by the limits x1, x2, y1, y2.
  Integration is performed by calling AdaptiveGaussianQuadrature recursively.*/

  long double DummyParam;
  FuncPtr = func;
  NbDigitsAGQ = digits;
  yy1 = c; yy2 = d;
  return AdaptiveGaussianQuadrature(AGQf1,a,b,NbDigitsAGQ,DummyParam);
}
 
double AGQf1(const double x){ 
  /*        /y2
    f1 = Int|  func(x,y)dy
            /y1        
  */

  long double DummyParam;
  xsav=x;
  return AdaptiveGaussianQuadrature(f2,yy1,yy2,NbDigitsAGQ,DummyParam);
}

double Error1D(const double x, const FunctionType1D func1, const FunctionType1D func2){
  /*returns the absolute difference between the evaluations of two different functions
  in the same point "x" of the 1D space*/

  return fabs(func1(x)-func2(x));
}


double error2D(const FunctionType2D func1, const FunctionType2D func2, const double x, const double y){
/*returns the absolute difference between the evaluations of two different functions
  in the same point "(x,y)" of the 2D space*/

  return fabs(func1(x,y)-func2(x,y)); 
}


/******************************************************************************************
Function for generating the geom coeff. for a cartesian cell
******************************************************************************************/

double GeomCoeffCartesian(int p1, int p2, double deltaX, double deltaY, double deltaXC, double deltaYC){

  /* p1 -> the first power coefficient
     p2 -> the second power coefficient
     deltaX -> the grid size in the X direction
     deltaY -> the grid size in the Y direction
     deltaXC -> the X distance between the center of the reconstructed cell and that of the cell used in the reconstruction
     deltaYC -> the Y distance between the center of the reconstructed cell and that of the cell used in the reconstruction

     Obs. To compute the coefficient of the reconstructed cell, deltaXC and deltaYC must be ZERO. 
 */

  double val1, val2;
  double coef_x1, coef_x2, coef_y1, coef_y2;

  val1 = val2 = 0.0;
  coef_x1 = deltaX/2  + deltaXC;
  coef_x2 = -deltaX/2 + deltaXC;
  coef_y1 = deltaY/2  + deltaYC;
  coef_y2 = -deltaY/2 + deltaYC;

  for (int m=1; m<=p1+1; m++){
    val1 += pow(coef_x1,p1+1-m)*pow(coef_x2,m-1);
  }
  for (int l=1; l<=p2+1; l++){
    val2 += pow(coef_y1,p2+1-l)*pow(coef_y2,l-1);
  }

  return val1*val2/((p1+1)*(p2+1));
}
