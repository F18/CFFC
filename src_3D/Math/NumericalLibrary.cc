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
 * \return the value of the integrals up to 4th-order (i.e. OrderX + OrderY <= 4)
 *
 *********************************************************************************/
double PolynomLineIntegration(const double & N1x, const double & N1y,
			      const double & N2x, const double & N2y,
			      const double & xCC, const double & yCC,
			      const int &OrderX,   const int &OrderY){
  
  double DX(N2x - N1x);
  double DY(N2y - N1y);

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
    case 4:
      /* OrderX=0, OrderY=4 */
      return ( DX * pow(DY, 0.5e1) / 0.6e1 + DY * (N1x * pow(DY, 0.4e1) - 0.4e1 * DX * pow(DY, 0.3e1) * yCC - 
						   xCC * pow(DY, 0.4e1) + 0.4e1 * DX * N1y * pow(DY, 0.3e1)) / 0.5e1 + 
	       DY * (-0.4e1 * xCC * N1y * pow(DY, 0.3e1) + 0.4e1 * xCC * pow(DY, 0.3e1) * yCC + 
		     0.6e1 * DX * N1y * N1y * DY * DY + 0.4e1 * N1x * N1y * pow(DY, 0.3e1) - 
		     0.12e2 * DX * N1y * DY * DY * yCC - 0.4e1 * N1x * pow(DY, 0.3e1) * yCC +
		     0.6e1 * DX * DY * DY * yCC * yCC) / 0.4e1 +
	       DY * (0.6e1 * N1x * DY * DY * yCC * yCC - 0.6e1 * xCC * N1y * N1y * DY * DY - 
		     0.6e1 * xCC * DY * DY * yCC * yCC + 0.4e1 * DX * pow(N1y, 0.3e1) * DY - 
		     0.4e1 * DX * DY * pow(yCC, 0.3e1) + 0.12e2 * DX * N1y * DY * yCC * yCC + 
		     0.6e1 * N1x * N1y * N1y * DY * DY - 0.12e2 * DX * N1y * N1y * DY * yCC - 
		     0.12e2 * N1x * N1y * DY * DY * yCC + 0.12e2 * xCC * N1y * DY * DY * yCC) / 0.3e1 + 
	       DY * (DX * pow(yCC, 0.4e1) - 0.4e1 * DX * N1y * pow(yCC, 0.3e1) + 0.6e1 * DX * N1y * N1y * yCC * yCC + 
		     0.4e1 * xCC * DY * pow(yCC, 0.3e1) - 0.12e2 * xCC * N1y * DY * yCC * yCC + 
		     0.12e2 * xCC * N1y * N1y * DY * yCC + 0.4e1 * N1x * pow(N1y, 0.3e1) * DY + 
		     DX * pow(N1y, 0.4e1) - 0.4e1 * xCC * pow(N1y, 0.3e1) * DY - 0.4e1 * DX * pow(N1y, 0.3e1) * yCC - 
		     0.4e1 * N1x * DY * pow(yCC, 0.3e1) - 0.12e2 * N1x * N1y * N1y * DY * yCC + 
		     0.12e2 * N1x * N1y * DY * yCC * yCC) / 0.2e1 + 
	       DY * (0.4e1 * xCC * N1y * pow(yCC, 0.3e1) + 0.6e1 * N1x * N1y * N1y * yCC * yCC + 
		     N1x * pow(yCC, 0.4e1) + 0.4e1 * xCC * pow(N1y, 0.3e1) * yCC - 
		     0.6e1 * xCC * N1y * N1y * yCC * yCC - 0.4e1 * N1x * pow(N1y, 0.3e1) * yCC - 
		     0.4e1 * N1x * N1y * pow(yCC, 0.3e1) - xCC * pow(yCC, 0.4e1) + N1x * pow(N1y, 0.4e1) - 
		     xCC * pow(N1y, 0.4e1)) );
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
    case 3:
      /* OrderX=1, OrderY=3 */
      return ((DX * DX / 0.6e1 + xCC * xCC / 0.4e1 - N1x * xCC / 0.2e1 + N1x * N1x / 0.4e1 - 
	       0.2e1 / 0.5e1 * DX * xCC + 0.2e1 / 0.5e1 * N1x * DX) * pow(DY, 0.4e1) + 
	      (-0.3e1 / 0.2e1 * N1x * DX * yCC + 0.3e1 / 0.2e1 * N1x * DX * N1y + 
	       0.3e1 / 0.2e1 * DX * xCC * yCC - 0.3e1 / 0.2e1 * DX * xCC * N1y + 
	       0.2e1 * N1x * xCC * yCC + N1x * N1x * N1y + xCC * xCC * N1y - N1x * N1x * yCC - 
	       xCC * xCC * yCC - 0.2e1 * N1x * xCC * N1y + 0.3e1 / 0.5e1 * DX * DX * N1y - 
	       0.3e1 / 0.5e1 * DX * DX * yCC) * pow(DY, 0.3e1) + 
	      (0.3e1 / 0.2e1 * N1x * N1x * N1y * N1y + 0.6e1 * N1x * xCC * N1y * yCC + 
	       0.3e1 / 0.2e1 * xCC * xCC * yCC * yCC - 0.3e1 * N1x * xCC * N1y * N1y + 
	       0.3e1 / 0.2e1 * xCC * xCC * N1y * N1y + 0.3e1 / 0.2e1 * N1x * N1x * yCC * yCC - 
	       0.3e1 * N1x * xCC * yCC * yCC - 0.3e1 * N1x * N1x * N1y * yCC - 0.3e1 * xCC * xCC * N1y * yCC + 
	       0.3e1 / 0.4e1 * DX * DX * yCC * yCC + 0.3e1 / 0.4e1 * DX * DX * N1y * N1y - 
	       0.3e1 / 0.2e1 * DX * DX * N1y * yCC + 0.4e1 * DX * xCC * N1y * yCC + 
	       0.2e1 * N1x * DX * N1y * N1y - 0.2e1 * DX * xCC * yCC * yCC + 0.2e1 * N1x * DX * yCC * yCC - 
	       0.2e1 * DX * xCC * N1y * N1y - 0.4e1 * N1x * DX * N1y * yCC) * DY * DY + 
	      (DX * xCC * pow(yCC, 0.3e1) - 0.3e1 * N1x * DX * N1y * N1y * yCC + 
	       0.3e1 * N1x * DX * N1y * yCC * yCC + 0.3e1 * DX * xCC * N1y * N1y * yCC + 
	       N1x * DX * pow(N1y, 0.3e1) - 0.3e1 * DX * xCC * N1y * yCC * yCC - 
	       N1x * DX * pow(yCC, 0.3e1) - DX * xCC * pow(N1y, 0.3e1) - DX * DX * pow(yCC, 0.3e1) / 0.3e1 - 
	       DX * DX * N1y * N1y * yCC + DX * DX * pow(N1y, 0.3e1) / 0.3e1 + DX * DX * N1y * yCC * yCC - 
	       xCC * xCC * pow(yCC, 0.3e1) + 0.3e1 * xCC * xCC * N1y * yCC * yCC - 
	       0.2e1 * N1x * xCC * pow(N1y, 0.3e1) + 0.6e1 * N1x * xCC * N1y * N1y * yCC - 
	       0.6e1 * N1x * xCC * N1y * yCC * yCC + 0.3e1 * N1x * N1x * N1y * yCC * yCC - 
	       0.3e1 * N1x * N1x * N1y * N1y * yCC - 0.3e1 * xCC * xCC * N1y * N1y * yCC + 
	       xCC * xCC * pow(N1y, 0.3e1) + 0.2e1 * N1x * xCC * pow(yCC, 0.3e1) + 
	       N1x * N1x * pow(N1y, 0.3e1) - N1x * N1x * pow(yCC, 0.3e1)) * DY );
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
    case 2:
      /* OrderX=2, OrderY=2 */
      return ((pow(DX, 0.3e1) / 0.6e1 + N1x * xCC * xCC - N1x * N1x * xCC - pow(xCC, 0.3e1) / 0.3e1 + 
	       pow(N1x, 0.3e1) / 0.3e1 - 0.3e1 / 0.2e1 * N1x * DX * xCC + 0.3e1 / 0.4e1 * N1x * N1x * DX + 
	       0.3e1 / 0.4e1 * DX * xCC * xCC + 0.3e1 / 0.5e1 * N1x * DX * DX - 
	       0.3e1 / 0.5e1 * DX * DX * xCC) * pow(DY, 0.3e1) +
	      (-0.3e1 / 0.2e1 * N1x * DX * DX * yCC - 0.3e1 / 0.2e1 * DX * DX * xCC * N1y + 0.3e1 / 0.2e1 * N1x * DX * DX * N1y + 
	       0.3e1 / 0.2e1 * DX * DX * xCC * yCC + 0.2e1 / 0.5e1 * pow(DX, 0.3e1) * N1y - 
	       0.2e1 / 0.5e1 * pow(DX, 0.3e1) * yCC - 0.3e1 * N1x * xCC * xCC * yCC - pow(N1x, 0.3e1) * yCC + 
	       0.3e1 * N1x * xCC * xCC * N1y + 0.3e1 * N1x * N1x * xCC * yCC + pow(xCC, 0.3e1) * yCC + 
	       pow(N1x, 0.3e1) * N1y - pow(xCC, 0.3e1) * N1y - 0.3e1 * N1x * N1x * xCC * N1y - 
	       0.2e1 * N1x * N1x * DX * yCC + 0.2e1 * N1x * N1x * DX * N1y + 0.2e1 * DX * xCC * xCC * N1y - 
	       0.4e1 * N1x * DX * xCC * N1y + 0.4e1 * N1x * DX * xCC * yCC - 0.2e1 * DX * xCC * xCC * yCC) * DY * DY +
	      (0.6e1 * N1x * N1x * xCC * N1y * yCC - 0.6e1 * N1x * xCC * xCC * N1y * yCC + 
	       0.3e1 / 0.2e1 * DX * xCC * xCC * N1y * N1y - pow(xCC, 0.3e1) * N1y * N1y + 
	       0.2e1 * DX * DX * xCC * N1y * yCC - 0.3e1 * DX * xCC * xCC * N1y * yCC - 
	       0.3e1 * N1x * N1x * xCC * N1y * N1y + 0.2e1 * pow(xCC, 0.3e1) * N1y * yCC + 
	       0.3e1 * N1x * xCC * xCC * yCC * yCC + 0.3e1 / 0.2e1 * N1x * N1x * DX * yCC * yCC - 
	       0.3e1 * N1x * N1x * xCC * yCC * yCC + 0.3e1 / 0.2e1 * N1x * N1x * DX * N1y * N1y + 
	       0.3e1 / 0.2e1 * DX * xCC * xCC * yCC * yCC - 0.2e1 * N1x * DX * DX * N1y * yCC + 
	       0.6e1 * N1x * DX * xCC * N1y * yCC - 0.3e1 * N1x * DX * xCC * yCC * yCC - 
	       0.3e1 * N1x * DX * xCC * N1y * N1y + N1x * DX * DX * yCC * yCC - DX * DX * xCC * N1y * N1y - 
	       0.2e1 * pow(N1x, 0.3e1) * N1y * yCC + 0.3e1 * N1x * xCC * xCC * N1y * N1y - 
	       DX * DX * xCC * yCC * yCC - pow(DX, 0.3e1) * N1y * yCC / 0.2e1 + pow(N1x, 0.3e1) * N1y * N1y - 
	       0.3e1 * N1x * N1x * DX * N1y * yCC + N1x * DX * DX * N1y * N1y + pow(N1x, 0.3e1) * yCC * yCC - 
	       pow(xCC, 0.3e1) * yCC * yCC + pow(DX, 0.3e1) * yCC * yCC / 0.4e1 + 
	       pow(DX, 0.3e1) * N1y * N1y / 0.4e1) * DY );
    } /* endswitch (OrderY) */
    break;

  case 3:
    switch(OrderY){
    case 0:
      /* OrderX=3, OrderY=0 */
      return DY*(DX*DX*DX*( 0.2*DX + N1x - xCC) + 
		 (0.2e1*DX*DX*xCC*xCC - 0.4e1*N1x*DX*DX*xCC + 0.2e1*N1x*N1x*DX*DX) + 
		 (-0.2e1*DX*xCC*xCC*xCC - 0.6e1*N1x*N1x*DX*xCC + 0.2e1*N1x*N1x*N1x*DX + 0.6e1*N1x*DX*xCC*xCC) + 
		 (N1x*N1x*N1x*(N1x - 0.4e1*xCC) + xCC*xCC*xCC*(xCC - 0.4e1*N1x)  + 0.6e1*N1x*N1x*xCC*xCC));
    case 1:
      /* OrderX=3, OrderY=1 */
      return ( (pow(DX, 0.4e1) / 0.6e1 - 0.4e1 / 0.5e1 * pow(DX, 0.3e1) * xCC + 0.4e1 / 0.5e1 * N1x * pow(DX, 0.3e1) +
		0.3e1 / 0.2e1 * N1x * N1x * DX * DX - 0.3e1 * N1x * DX * DX * xCC + 0.3e1 / 0.2e1 * DX * DX * xCC * xCC +
		0.4e1 / 0.3e1 * pow(N1x, 0.3e1) * DX - 0.4e1 / 0.3e1 * DX * pow(xCC, 0.3e1) + 
		0.4e1 * N1x * DX * xCC * xCC - 0.4e1 * N1x * N1x * DX * xCC + pow(xCC, 0.4e1) / 0.2e1 - 
		0.2e1 * pow(N1x, 0.3e1) * xCC + 0.3e1 * N1x * N1x * xCC * xCC - 0.2e1 * N1x * pow(xCC, 0.3e1) + 
		pow(N1x, 0.4e1) / 0.2e1) * DY * DY + 
	       (0.2e1 * N1x * N1x * DX * DX * N1y - 0.4e1 * N1x * pow(xCC, 0.3e1) * N1y + 
		pow(N1x, 0.4e1) * N1y - pow(N1x, 0.4e1) * yCC + 0.4e1 * pow(N1x, 0.3e1) * xCC * yCC - 
		0.2e1 * DX * pow(xCC, 0.3e1) * N1y + pow(xCC, 0.4e1) * N1y - 
		pow(xCC, 0.4e1) * yCC + pow(DX, 0.4e1) * N1y / 0.5e1 - 0.4e1 * N1x * DX * DX * xCC * N1y + 
		0.6e1 * N1x * DX * xCC * xCC * N1y + 0.4e1 * N1x * pow(xCC, 0.3e1) * yCC + 
		0.2e1 * DX * DX * xCC * xCC * N1y + 0.2e1 * pow(N1x, 0.3e1) * DX * N1y - 
		0.4e1 * pow(N1x, 0.3e1) * xCC * N1y - 0.6e1 * N1x * DX * xCC * xCC * yCC - 
		0.6e1 * N1x * N1x * DX * xCC * N1y + 0.6e1 * N1x * N1x * DX * xCC * yCC - 
		0.2e1 * N1x * N1x * DX * DX * yCC - 0.2e1 * DX * DX * xCC * xCC * yCC + 
		0.2e1 * DX * pow(xCC, 0.3e1) * yCC - pow(DX, 0.3e1) * xCC * N1y + pow(DX, 0.3e1) * xCC * yCC - 
		pow(DX, 0.4e1) * yCC / 0.5e1 + N1x * pow(DX, 0.3e1) * N1y - N1x * pow(DX, 0.3e1) * yCC + 
		0.4e1 * N1x * DX * DX * xCC * yCC + 0.6e1 * N1x * N1x * xCC * xCC * N1y - 
		0.6e1 * N1x * N1x * xCC * xCC * yCC - 0.2e1 * pow(N1x, 0.3e1) * DX * yCC) * DY );
    }
    break;

  case 4:
    if (OrderY == 0){
      /* OrderX=4, OrderY=0 */
      return (pow(DX, 0.5e1) / 0.6e1 + N1x * pow(DX, 0.4e1) - pow(DX, 0.4e1) * xCC + 
	      0.5e1 / 0.2e1 * N1x * N1x * pow(DX, 0.3e1) - 0.5e1 * N1x * pow(DX, 0.3e1) * xCC + 
	      0.5e1 / 0.2e1 * pow(DX, 0.3e1) * xCC * xCC + 0.10e2 / 0.3e1 * pow(N1x, 0.3e1) * DX * DX - 
	      0.10e2 * N1x * N1x * DX * DX * xCC + 0.10e2 * N1x * DX * DX * xCC * xCC - 
	      0.10e2 / 0.3e1 * DX * DX * pow(xCC, 0.3e1) + 0.5e1 / 0.2e1 * pow(N1x, 0.4e1) * DX - 
	      0.10e2 * N1x * DX * pow(xCC, 0.3e1) - 0.10e2 * pow(N1x, 0.3e1) * DX * xCC + 
	      0.5e1 / 0.2e1 * DX * pow(xCC, 0.4e1) + 0.15e2 * N1x * N1x * DX * xCC * xCC + 
	      0.5e1 * N1x * pow(xCC, 0.4e1) + pow(N1x, 0.5e1) + 0.10e2 * pow(N1x, 0.3e1) * xCC * xCC - 
	      0.10e2 * N1x * N1x * pow(xCC, 0.3e1) - 0.5e1 * pow(N1x, 0.4e1) * xCC - pow(xCC, 0.5e1)) * DY;
    }
  } /* endswitch(OrderX) */

  std::cout << "PolynomLineIntegration ERROR: Power higher than the maximum allowed!\n";
  return 0.0;

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
 * \return the value of the integrals up to 4th-order (i.e. OrderX + OrderY <= 4)
 *
 * This subroutine is different than the PolynomLineIntegration() one in that
 * it tries to be more accurate for line segments characterized by large ratios
 * between the coordinates of the end points relative to the size of the segment.
 * It does that by exploiting the fact that the difference between the point coordinates
 * and the centroid might have the same order of magnitude as the segment sizes DX and DY.
 * Thus, round off errors are avoided.
 *********************************************************************************/
double PolynomLineIntegration2(const double & N1x, const double & N1y,
			       const double & N2x, const double & N2y,
			       const double & xCC, const double & yCC,
			       const int &OrderX,   const int &OrderY){
  
  double DX(N2x - N1x);		// DeltaX of the segment line
  double DY(N2y - N1y);		// DeltaY of the segment line
  double dXc(N1x - xCC);	// x-distance to centroid
  double dYc(N1y - yCC);	// y-distance to centroid

  // Temporary variables;
  double t1, t2, t3, t4, t5, t6, t7, t8;

  switch(OrderX){

  case 0:
    switch(OrderY){
    case 0:
      /* OrderX=0, OrderY=0 */
      return (0.5*DX + dXc)*DY;

    case 1:
      /* OrderX=0, OrderY=1 */
      t1 = DY*DY;
      return (DX*ONETHIRD + 0.5*dXc)*t1+(0.5*DX + dXc)*dYc*DY; 
  
    case 2:
      /* OrderX=0, OrderY=2 */
      t1 = DY*DY;
      t2 = dYc*dYc;
      return (dXc*ONETHIRD+0.25*DX)*DY*t1+(dXc+TWOTHIRDS*DX)*dYc*t1+(dXc+DX*0.5)*t2*DY; 

    case 3:
      /* OrderX=0, OrderY=3 */
      t1 = DY*DY;
      t2 = t1*t1;
      t3 = dYc*dYc;
      t4 = dYc*t3;
      return (DX*0.2+dXc*0.25)*t2+(dXc+0.75*DX)*dYc*DY*t1+(1.5*dXc+DX)*t3*t1+(DX*0.5+dXc)*t4*DY; 

    case 4:
      /* OrderX=0, OrderY=4 */
      t1 = DY*DY;
      t2 = t1*t1;
      t3 = dYc*dYc;
      t4 = dYc*t3;
      t5 = t3*t3;
      return ( (dXc*0.2+DX*ONESIXTH)*DY*t2 + (dXc+0.8*DX)*dYc*t2 + (2.0*dXc+1.5*DX)*t3*DY*t1 +
	       (2.0*dXc+FOURTHIRDS*DX)*t4*t1+(dXc+DX*0.5)*t5*DY );
    } /* endswitch (OrderY) */
    break;

  case 1:
    switch(OrderY){
    case 0:
      /* OrderX=1, OrderY=0 */
      t1 = DX*DX;
      t2 = dXc*dXc;
      return (t1*ONETHIRD+dXc*DX+t2)*DY;
  
    case 1:
      /* OrderX=1, OrderY=1 */
      t1 = DX*DX;
      t2 = dXc*DX;
      t3 = dXc*dXc;
      t4 = DY*DY;
      return (t1*0.25+TWOTHIRDS*t2+t3*0.5)*t4 + (t1*ONETHIRD+t2+t3)*dYc*DY;
      
    case 2:
      /* OrderX=1, OrderY=2 */
      t1 = dXc*dXc;
      t2 = dXc*DX;
      t3 = DX*DX;
      t4 = DY*DY;
      t5 = dYc*dYc;
      return ((t1*ONETHIRD+t2*0.5+t3*0.2)*DY + (t3*0.5+FOURTHIRDS*t2+t1)*dYc)*t4+(t1+t2+t3*ONETHIRD)*t5*DY;

    case 3:
      /* OrderX=1, OrderY=3 */
      t1 = DX*DX;
      t2 = dXc*dXc;
      t3 = dXc*DX;
      t4 = DY*DY;
      t5 = t4*t4;
      t6 = dYc*dYc;
      t7 = dYc*t6;
      return ( (t1*ONESIXTH+t2*0.25+0.4*t3)*t5 + (1.5*t3+t2+0.6*t1)*dYc*DY*t4+
	       (1.5*t2+0.75*t1+2.0*t3)*t6*t4 + (t3+t1*ONETHIRD+t2)*t7*DY );
    } /* endswitch (OrderY) */
    break;

  case 2:
    switch(OrderY){
    case 0:
      /* OrderX=2, OrderY=0 */
      t1 = DX*DX;
      t2 = dXc*dXc;
      return ((DX*0.25+dXc)*t1 + (1.5*DX+dXc)*t2)*DY;
  
    case 1:
      /* OrderX=2, OrderY=1 */
      t1 = DX*DX;
      t2 = DX*t1;
      t3 = dXc*t1;
      t4 = dXc*dXc;
      t5 = t4*DX;
      t6 = dXc*t4;
      t7 = DY*DY;
      return (t2*0.2+0.75*t3+t5+t6*0.5)*t7+(t2*0.25+t3+1.5*t5+t6)*dYc*DY;

    case 2:
      /* OrderX=2, OrderY=2 */
      t1 = DX*DX;
      t2 = DX*t1;
      t3 = dXc*dXc;
      t4 = dXc*t3;
      t5 = t3*DX;
      t6 = dXc*t1;
      t7 = DY*DY;
      t8 = dYc*dYc;
      return ( ( (t2*ONESIXTH+t4*ONETHIRD+0.75*t5+0.6*t6)*DY + (1.5*t6+0.4*t2+t4+2.0*t5)*dYc)*t7+
	       (1.5*t5+t6+t2*0.25+t4)*t8*DY );
    } /* endswitch (OrderY) */
    break;

  case 3:
    switch(OrderY){
    case 0:
      /* OrderX=3, OrderY=0 */
      t1 = DX*DX;
      t2 = t1*t1;
      t3 = dXc*dXc;
      t4 = t3*t3;
      return (t2*0.2 + dXc*DX*t1 + 2.0*t3*(t1+dXc*DX) +t4)*DY;

    case 1:
      /* OrderX=3, OrderY=1 */
      t1 = DX*DX;
      t2 = t1*t1;
      t3 = dXc*DX*t1;
      t4 = dXc*dXc;
      t5 = t4*t1;
      t6 = dXc*t4*DX;
      t7 = t4*t4;
      t8 = DY*DY;
      return ( (t2*ONESIXTH+0.8*t3+1.5*t5+FOURTHIRDS*t6+t7*0.5)*t8 + (t2*0.2+t3+2.0*t5+2.0*t6+t7)*dYc*DY );
    }
    break;

  case 4:
    if (OrderY == 0){
      /* OrderX=4, OrderY=0 */
      t1 = DX*DX;
      t2 = t1*t1;
      t3 = dXc*dXc;
      t4 = t3*t3;
      return (DX*(t2*ONESIXTH+2.5*(t3*t1+t4)) + dXc*(t2+TENTHIRDS*t3*t1+t4) )*DY;
    }
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
  //int TRUE(1);

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
    odd=true;
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

//! Abscissae for 3-point Gaussian (i.e. 0.5*(1 - sqrt(3)/sqrt(5)), 0.5*(1-0) , 0.5*(1 + sqrt(3)/sqrt(5)) )
const double GaussQuadratureData::GQ3_Abscissa[3] = {0.1127016653792583114820736,
						     0.5 ,
						     0.8872983346207416885179264};
//! Weights for 3-point Gaussian (i.e. 0.5*5/9, 0.5*8/9, 0.5*5/9)
const double GaussQuadratureData::GQ3_Weight[3] = {2.77777777777778e-1,
						   4.44444444444444e-1,
						   2.77777777777778e-1};

//! Abscissae for 5-point Gaussian (i.e. 0.5*(1 - 0.90617984593866399280),
//                                       0.5*(1 - 0.53846931010568309104),
//                                       0.5*(1 - 0),
//                                       0.5*(1 + 0.53846931010568309104),
//                                       0.5*(1 + 0.90617984593866399280) )
const double GaussQuadratureData::GQ5_Abscissa[5] = {4.6910077030668e-2,
						     2.30765344947158e-1,
						     0.5,
						     7.69234655052842e-1,
						     9.53089922969332e-1};
//! Weights for 5-point Gaussian (i.e. 0.5*0.23692688505618908751,
//                                     0.5*0.47862867049936646804,
//                                     0.5*0.56888888888888888889,
//                                     0.5*0.47862867049936646804,
//                                     0.5*0.23692688505618908751)
const double GaussQuadratureData::GQ5_Weight[5] = {1.18463442528094543755e-1,
						   2.39314335249683234020e-1,
						   2.84444444444444444445e-1,
						   2.39314335249683234020e-1,
						   1.18463442528094543755e-1};


