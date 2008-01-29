/*!\file NumericalLibrary.cc
  \brief Source file providing implementation of subroutines prototyped in NumericalLibrary.h */

/* Include required C++ libraries. */
#include <complex>
#include <limits>

/* Using std namespace functions */
using std::complex;

/* Include CFFC header files */
#include "NumericalLibrary.h"

/***********************************************************************//**
 * Compute the integral \f$ I = \int x dy \f$ along a segment line.
 * This integral is computed with an analytic expression.
 * This function arises in the area computation of plan figures 
 * with curved boundaries when the curves are approximated by segment lines.
 * For the double (area) integral, the powers of the general polynomial
 * function in "x" and "y" are ZERO.    
 *
 * \param N1x the x-coordinate of the segment first end point.
 * \param N1y the y-coordinate of the segment first end point.
 * \param N2x the x-coordinate of the segment second end point.
 * \param N2y the y-coordinate of the segment second end point.
 * \return the value of the integral
 ************************************************************************/
double ZeroLineIntegration(const double & N1x, const double & N1y,
			   const double & N2x, const double & N2y){

  return (N2y - N1y)*(N1x + 0.5*(N2x - N1x));
}

/***********************************************************************//**
 * Compute the integral \f$ I = \int (x - xc)^{(OrderX+1)} * (y - yc)^OrderY dy \f$
 * along a segment line.
 * The result of this polynomial function integration is determined with an analytic expression.
 *
 * \param N1x the x-coordinate of the line first end point
 * \param N1y the y-coordinate of the line first end point
 * \param N2x the x-coordinate of the line second end point
 * \param N2y the y-coordinate of the line second end point
 * \param xCC the xc-coordinate
 * \param yCC the yc-coordinate
 * \return the value of the integrals up to 3rd-order (i.e. OrderX + OrderY <= 3)
 *
 * \todo extend the procedure up to 4th-order (required for viscous terms)!
 *********************************************************************************/
double PolynomLineIntegration(const double & N1x, const double & N1y,
			      const double & N2x, const double & N2y,
			      const double & xCC, const double & yCC,
			      const int &OrderX,   const int &OrderY){
  
  double DX = N2x - N1x;
  double DY = N2y - N1y;

  switch(OrderX){

  case 0:
    switch(OrderY){
    case 0:
      /* OrderX=0, OrderY=0 */
      return DY*(N1x + DX*0.5 - xCC);

    case 1:
      /* OrderX=0, OrderY=1 */
      return DY*(DX*DY/0.3e1 + (N1y-yCC)*(N1x - xCC + 0.5*DX) + 0.5*DY*(N1x-xCC));
  
    case 2:
      /* OrderX=0, OrderY=2 */
      return DY*(DX*DY*DY*0.25 + DY*(0.2e1*DX*(N1y - yCC) + DY*(N1x - xCC))/0.3e1 + 
		 (0.2e1*DY*(xCC*yCC - xCC*N1y + N1x*N1y - N1x*yCC) + DX*(N1y-yCC)*(N1y-yCC) )*0.5 + 
		 (N1x-xCC)*(N1y*N1y + yCC*yCC) + 0.2e1*N1y*yCC*(xCC - N1x) );

    case 3:
      /* OrderX=0, OrderY=3 */
      return DY*(DX*DY*DY*DY*0.2 + (-0.3e1*DX*DY*DY*yCC - xCC*DY*DY*DY + N1x*DY*DY*DY + 0.3e1*DX*N1y*DY*DY)*0.25 + 
		 (-0.3e1*xCC*N1y*DY*DY + 0.3e1*DX*DY*yCC*yCC + 0.3e1*DX*N1y*N1y*DY + 0.3e1*N1x*N1y*DY*DY +
		  0.3e1*xCC*DY*DY*yCC - 0.3e1*N1x*DY*DY*yCC - 0.6e1*DX*N1y*DY*yCC)/0.3e1 + 
		 (-0.3e1*xCC*DY*yCC*yCC + 0.3e1*DX*N1y*yCC*yCC + 0.3e1*N1x*N1y*N1y*DY + DX*N1y*N1y*N1y + 
		  0.3e1*N1x*DY*yCC*yCC - DX*yCC*yCC*yCC + 0.6e1*xCC*N1y*DY*yCC - 0.3e1*DX*N1y*N1y*yCC - 
		  0.3e1*xCC*N1y*N1y*DY - 0.6e1*N1x*N1y*DY*yCC)*0.5 + 
		 (N1x*N1y*N1y*N1y - N1x*yCC*yCC*yCC - xCC*N1y*N1y*N1y + xCC*yCC*yCC*yCC + 0.3e1*xCC*N1y*N1y*yCC - 
		  0.3e1*xCC*N1y*yCC*yCC - 0.3e1*N1x*N1y*N1y*yCC + 0.3e1*N1x*N1y*yCC*yCC));
    } /* endswitch (OrderY) */
    break;

  case 1:
    switch(OrderY){
    case 0:
      /* OrderX=1, OrderY=0 */
      return DY*(DX*DX/0.3e1 -DX*xCC + N1x*DX + N1x*N1x + xCC*xCC - 0.2e1*N1x*xCC);
  
    case 1:
      /* OrderX=1, OrderY=1 */
      return DY*(DX*DX*DY*0.25 + (DX*DX*N1y - 0.2e1*DX*xCC*DY - DX*DX*yCC + 0.2e1*N1x*DX*DY)/0.3e1 + 
		 (N1x*N1x*DY + xCC*xCC*DY)*0.5 - N1x*DX*yCC + N1x*DX*N1y + DX*xCC*yCC - N1x*xCC*DY - DX*xCC*N1y + 
		 (N1x*N1x*N1y - N1x*N1x*yCC + xCC*xCC*N1y - xCC*xCC*yCC + 0.2e1*N1x*xCC*(yCC- N1y)));
      
    case 2:
      /* OrderX=1, OrderY=2 */
      return DY*(DX*DX*DY*DY*0.2 + 
		 (N1x*DX*DY*DY - DX*xCC*DY*DY - DX*DX*DY*yCC + DX*DX*N1y*DY)*0.5 + 
		 (-0.2e1*N1x*xCC*DY*DY + N1x*N1x*DY*DY + 0.4e1*N1x*DX*N1y*DY + DX*DX*N1y*N1y - 0.4e1*DX*xCC*N1y*DY - 
		  0.2e1*DX*DX*N1y*yCC + xCC*xCC*DY*DY - 0.4e1*N1x*DX*DY*yCC + 0.4e1*DX*xCC*DY*yCC + DX*DX*yCC*yCC)/0.3e1 +
		 (-DX*xCC*N1y*N1y - 0.2e1*N1x*xCC*N1y*DY + 0.2e1*N1x*xCC*DY*yCC + N1x*DX*yCC*yCC + 
		  N1x*N1x*N1y*DY - DX*xCC*yCC*yCC - xCC*xCC*DY*yCC + xCC*xCC*N1y*DY + 
		  N1x*DX*N1y*N1y - 0.2e1*N1x*DX*N1y*yCC + 0.2e1*DX*xCC*N1y*yCC - N1x*N1x*DY*yCC) + 
		 (N1x*N1x*yCC*yCC + xCC*xCC*N1y*N1y + xCC*xCC*yCC*yCC - 0.2e1*N1x*xCC*N1y*N1y - 0.2e1*N1x*xCC*yCC*yCC + 
		  N1x*N1x*N1y*N1y - 0.2e1*N1x*N1x*N1y*yCC - 0.2e1*xCC*xCC*N1y*yCC + 0.4e1*N1x*xCC*N1y*yCC));
    } /* endswitch (OrderY) */
    break;

  case 2:
    switch(OrderY){
    case 0:
      /* OrderX=2, OrderY=0 */
      return DY*(DX*DX*DX*0.25 + DX*DX*(N1x - xCC) + 0.15e1*DX*(-0.2e1*N1x*xCC + xCC*xCC + N1x*N1x) + N1x*N1x*N1x +
		 0.3e1*N1x*xCC*(xCC - N1x) - xCC*xCC*xCC );
  
    case 1:
      /* OrderX=2, OrderY=1 */
      return DY*(DX*DX*DX*DY*0.2 + (0.3e1*DX*DX*DY*(N1x - xCC) + DX*DX*DX*(N1y - yCC))*0.25 + 
		 (N1x*DX*DX*(N1y - yCC) + DX*DX*xCC*(yCC - N1y) + DX*xCC*DY*(xCC - 0.2e1*N1x) + 
		  N1x*N1x*DX*DY) + 
		 (N1x*N1x*N1x*DY - xCC*xCC*xCC*DY - 0.6e1*N1x*DX*xCC*N1y + 0.3e1*N1x*xCC*xCC*DY - 
		  0.3e1*N1x*N1x*DX*yCC - 0.3e1*N1x*N1x*xCC*DY - 0.3e1*DX*xCC*xCC*yCC + 0.3e1*DX*xCC*xCC*N1y + 
		  0.6e1*N1x*DX*xCC*yCC + 0.3e1*N1x*N1x*DX*N1y)*0.5 + 
		 (xCC*xCC*xCC*(yCC - N1y) + 0.3e1*N1x*N1x*xCC*(yCC - N1y) - 0.3e1*N1x*xCC*xCC*(yCC - N1y)
		  + N1x*N1x*N1x*(N1y - yCC)));
  
    } /* endswitch (OrderY) */
    break;

  case 3:
    if (OrderY == 0){
      /* OrderX=3, OrderY=0 */
      return DY*(DX*DX*DX*( 0.2*DX + N1x - xCC) + 
		 (0.2e1*DX*DX*xCC*xCC - 0.4e1*N1x*DX*DX*xCC + 0.2e1*N1x*N1x*DX*DX) + 
		 (-0.2e1*DX*xCC*xCC*xCC - 0.6e1*N1x*N1x*DX*xCC + 0.2e1*N1x*N1x*N1x*DX + 0.6e1*N1x*DX*xCC*xCC) + 
		 (N1x*N1x*N1x*(N1x - 0.4e1*xCC) + xCC*xCC*xCC*(xCC - 0.4e1*N1x)  + 0.6e1*N1x*N1x*xCC*xCC));
    }
    break;
  } /* endswitch(OrderX) */

  std::cout << "PolynomLineIntegration ERROR: Power higher than the maximum allowed!\n";
  return 0.0;

}

// frenel()
void frenel(double x, double &s, double &c)
{
  /* returns the sin (s) and cos (c) values for the Frenel function applied to "x"*/
  double EPS(6.0e-8);
  int MAXIT(100);
  double XMIN(1.5), PIBY2(PI/2.0), FPMIN(1.0e-30);
  int TRUE(1);

  int k,n,odd;
  float a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
  complex<double> b, cc, d, h, del, cs, a_compl;
  
  ax=fabs(x);
  if (ax < sqrt(FPMIN)) {
    s=0.0;
    c=ax;
  } else if (ax <= XMIN) {
    sum=sums=0.0;
    sumc=ax;
    sign=1.0;
    fact=PIBY2*ax*ax;
    odd=TRUE;
    term=ax;
    n=3;
    for (k=1; k<=MAXIT; ++k) {
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
    s=sums;
    c=sumc;
  } else {
    pix2 = PI*ax*ax;
    b = complex<double>(1.0,-pix2);
    cc = 1.0/FPMIN;
    d = h = 1.0/b;
    n=-1;
    for (k=2; k<=MAXIT; ++k){
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
    cs = ( complex<double>(0.5,0.5)*(complex<double>(1.0,0.0)-h*complex<double>(cos(0.5*pix2),sin(0.5*pix2))) );
    c = cs.real();
    s = cs.imag();
  }
  if (x < 0.0) {
    c = -c;
    s = -s;
  }
  return;
}

/* qgauss10() :returns the integral of the function "func" between a and b,
   by ten-point Gauss-Legendre integration: the function is evaluated exactly
   ten times at interior points in the range of integration.
   The integral is exact for polynomials up to order of 19. 
   Implementation based on the subroutine from Numerical Recipes*/
double qgauss10(const FunctionType1D func, const double a, const double b){

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

// Static variables:
static double xsav,yy1,yy2;
static FunctionType2D FuncPtr;
static int NbDigitsAGQ;

/* quad2d() : returns the integral of a user-supplied func over a two-dimensional 
   rectangular region, specified by the limits x1, x2, y1, y2.
   Integration is performed by calling qgauss10 recursively.*/
double quad2d(const FunctionType2D func, const double a, const double b, const double c, const double d){

  FuncPtr = func;
  yy1 = c; yy2 = d;
  return qgauss10(f1,a,b);
}
 
double f1(const double x){ 
  /*        /y2
	    f1 = Int|  func(x,y)dy
            /y1        
  */

  xsav=x;
  return qgauss10(f2,yy1,yy2);
}

double f2(const double y){
  // The integrand f(x,y) evaluated at fixed x and y
  return FuncPtr(xsav,y);
}

double quad2dAdaptiveGaussianQuadrature(const FunctionType2D func, const double a,
					const double b, const double c, const double d, int digits){
  /*returns the integral of a user-supplied func over a two-dimensional 
    rectangular region, specified by the limits x1, x2, y1, y2.
    Integration is performed by calling AdaptiveGaussianQuadrature recursively.*/

  long double DummyParam;
  FuncPtr = func;
  NbDigitsAGQ = digits;
  yy1 = c; yy2 = d;
  return AdaptiveGaussianQuadrature(AGQf1,a,b,DummyParam,NbDigitsAGQ);
}
 
double AGQf1(const double x){ 
  /*        /y2
	    f1 = Int|  func(x,y)dy
            /y1        
  */

  long double DummyParam;
  xsav=x;
  return AdaptiveGaussianQuadrature(f2,yy1,yy2,DummyParam,NbDigitsAGQ);
}

// ===  Static member variables GaussQuadratureData class ===
const double GaussQuadratureData::GQ1_Abscissa[1] = {0.5};	//!< Abscissa for 1-point Gaussian
const double GaussQuadratureData::GQ1_Weight[1] = {1.0};      //!< Weight for 1-point Gaussian

//! Abscissae for 2-point Gaussian (i.e. 0.5*(1 - 1/sqrt(3)), 0.5*(1 + 1/sqrt(3)))
const double GaussQuadratureData::GQ2_Abscissa[2] = {0.2113248654051871177454256,
						     0.7886751345948128822545744};
//! Weights for 2-point Gaussian
const double GaussQuadratureData::GQ2_Weight[2] = {0.5,
						   0.5};

//! Abscissae for 3-point Gaussian (i.e. 0.5*(1 - sqrt(3)/sqrt(5)), 0 , 0.5*(1 + 1/sqrt(3)) )
const double GaussQuadratureData::GQ3_Abscissa[3] = {0.1127016653792583114820736,
						     0.5 ,
						     0.8872983346207416885179264};
//! Weights for 3-point Gaussian (i.e. 0.5*5/9, 0.5*8/9, 0.5*5/9)
const double GaussQuadratureData::GQ3_Weight[3] = {2.77777777777778e-1,
						   4.44444444444444e-1,
						   2.77777777777778e-1};

