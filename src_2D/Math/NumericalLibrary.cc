/*!\file NumericalLibrary.cc
  \brief Source file providing implementation of subroutines prototyped in NumericalLibrary.h */

/* Include required C++ libraries. */
#include <complex>
#include <limits>

/* Using std namespace functions */
using std::complex;

/* Include CFDkit+caboodle header files */
#include "NumericalLibrary.h"

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
  long double LengthInterval = fabs(b - a); 
  long double StartPoint = a;
  long double EndPoint = b;
  long double SubDivisionPoint;
  double Division = 0.5;

  try{
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

  return 0.0;			// only for the compiler
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

//MakeReconstructionStencil(int,int,vector<int>) function
// Set the stencil for the 1D DD_ENO reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, vector<int> & i_index){

  // Obs. The first position (i_index[0]) corresponds to (iCell)

  switch(rings){

  case 3: // three rings of cells around (iCell)
    /* Third ring */
    i_index[5]=iCell-3;
    i_index[6]=iCell+3;

  case 2: // two rings of cells around (iCell)
    /* Second ring */
    i_index[3]=iCell-2;
    i_index[4]=iCell+2;

  case 1: // one ring of cells around (iCell,jCell)
    /* First ring */
    i_index[1]=iCell-1;
    i_index[2]=iCell+1;

  case 0: 
    i_index[0]=iCell; /* cell (iCell) */

  default: // general expression
    i_index[0]=iCell; /* cell (iCell) */
    for (int i=iCell-rings, Pos=1; i<=iCell+rings; ++i){
      if(i!=iCell){
	i_index[Pos] = i;
	++Pos;
      }
    }
  }
}

//MakeReconstructionStencil(int,int,int,vector<int>,vector<int>) function
// Set the stencil for the 2D DD_ENO reconstruction
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       vector<int> & i_index, vector<int> & j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    i_index[9] =iCell-2;  j_index[9]=jCell-2;
    i_index[10]=iCell-1; j_index[10]=jCell-2;
    i_index[11]=iCell  ; j_index[11]=jCell-2;
    i_index[12]=iCell+1; j_index[12]=jCell-2;
    i_index[13]=iCell+2; j_index[13]=jCell-2;
    i_index[14]=iCell-2; j_index[14]=jCell-1;
    i_index[15]=iCell+2; j_index[15]=jCell-1;
    i_index[16]=iCell-2; j_index[16]=jCell;
    i_index[17]=iCell+2; j_index[17]=jCell;
    i_index[18]=iCell-2; j_index[18]=jCell+1;
    i_index[19]=iCell+2; j_index[19]=jCell+1;
    i_index[20]=iCell-2; j_index[20]=jCell+2;
    i_index[21]=iCell-1; j_index[21]=jCell+2;
    i_index[22]=iCell  ; j_index[22]=jCell+2;
    i_index[23]=iCell+1; j_index[23]=jCell+2;
    i_index[24]=iCell+2; j_index[24]=jCell+2;

  case 1: // one ring of cells around (iCell,jCell)

    i_index[0]=iCell;   j_index[0]=jCell; /* cell (iCell,jCell) */
    /* First ring */
    i_index[1]=iCell-1; j_index[1]=jCell-1;
    i_index[2]=iCell;   j_index[2]=jCell-1;
    i_index[3]=iCell+1; j_index[3]=jCell-1;
    i_index[4]=iCell-1; j_index[4]=jCell;
    i_index[5]=iCell+1; j_index[5]=jCell;
    i_index[6]=iCell-1; j_index[6]=jCell+1;
    i_index[7]=iCell;   j_index[7]=jCell+1;
    i_index[8]=iCell+1; j_index[8]=jCell+1;
    break;

  default: // general expression
    i_index[0] = iCell;
    j_index[0] = jCell;
    for (int i=iCell-rings, Poz=1; i<=iCell+rings; ++i)
      for (int j=jCell-rings; j<=jCell+rings; ++j){
	if(!((i==iCell)&&(j==jCell)) ){
	  i_index[Poz] = i;
	  j_index[Poz] = j;
	  ++Poz;
	}
      }
  }//endswitch
 
}

//MakeReconstructionStencil(int,int,int,int,,vector<int>,vector<int>) function
// Enlarge the stencil for the 2D DD_ENO reconstruction used at curved boundaries
void MakeReconstructionStencil(const int & rings, const int & iCell, const int & jCell,
			       const int NorthCurvedBnd, const int SouthCurvedBnd,
			       const int EastCurvedBnd, const int WestCurvedBnd,
			       const int &ICl, const int &ICu, const int &JCl, const int &JCu,
			       int & StencilDimension, 
			       vector<int> & i_index, vector<int> & j_index){

  // Obs. The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell)

  i_index[0] = iCell;
  j_index[0] = jCell;
  int i,j, Imin, Imax, Jmin, Jmax, Poz;

  /* Determine Imin, Imax, Jmin, Jmax */
  Imin = iCell-rings; Imax = iCell+rings;
  Jmin = jCell-rings; Jmax = jCell+rings;

  if( NorthCurvedBnd==1 && Jmax > JCu ){ /* the north boundary is curved */
    if (jCell == JCu){
      Jmin -= 1;		/* add one extra layer in the j direction*/
      Jmax = JCu;
    } else {
      Jmax = JCu;		/* limit Jmax */
    }
  }

  if(SouthCurvedBnd==1 && Jmin < JCl){ 	/* the south boundary is curved */
    if (jCell == JCl){
      Jmax += 1;		/* add one extra layer in the j direction*/
      Jmin = JCl;
    } else {
      Jmin = JCl;		/* limit Jmin */
    }
  }

  if(EastCurvedBnd==1 && Imax > ICu){  /* the east boundary is curved */
    if (iCell == ICu){
      Imin -= 1;                /* add one extra layer in the i direction*/
      Imax = ICu;
    } else {
      Imax = ICu;		/* limit Imax */
    }
  }

  if(WestCurvedBnd==1 && Imin < ICl){  	/* the west boundary is curved */
    if (iCell == ICl){
      Imax += 1;                /* add one extra layer in the i direction*/
      Imin = ICl;
    } else {
      Imin = ICl;		/* limit Imin */
    }
  }

  /* Form stencil */
  for (i=Imin, Poz=1; i<=Imax; ++i)
    for (j=Jmin; j<=Jmax; ++j){
      if(!((i==iCell)&&(j==jCell)) ){
	i_index[Poz] = i;
	j_index[Poz] = j;
	++Poz;
      }//endif
    }// endfor

  StencilDimension = Poz;
}

/*************************************************************************
 * ZeroLineIntegration(StartPoint,EndPoint)                              *
 * Computes the analytical integral of the function with the form (x dy) *
 * along a segment line. This function arises in the computation of area *
 * of plan figures with curved boundaries when the curves are            *
 * approximated by segment lines. For the double (area) integral, the    *
 * powers of the general polynomial function in "x" and "y" are ZERO.    *
 * Input data:                                                           *
 *    - StartPoint.X (N1x)                                               *
 *    - StartPoint.Y (N1y)                                               *
 *    - EndPoint.X (N2x)                                                 *
 *    - EndPoint.Y (N2y)                                                 *
 * Returned data: -the value of the integral                             *
 ************************************************************************/
double ZeroLineIntegration(const double & N1x, const double & N1y,
			   const double & N2x, const double & N2y){

  return (N2y - N1y)*(N1x + 0.5*(N2x - N1x));
}

/*************************************************************************
 * PolynomLineIntegration(StartPoint,EndPoint,Centroid)                  *
 * Computes the analytical integral of polynomial functions with the     *
 * form (x - Centroid.X)^(OrderX+1) * (y - Centroid.Y)^OrderY along a    *
 * segment line.                                                         *
 * Input data:                                                           *
 *    - StartPoint.X (N1x)                                               *
 *    - StartPoint.Y (N1y)                                               *
 *    - EndPoint.X (N2x)                                                 *
 *    - EndPoint.Y (N2y)                                                 *
 *    - Centroid.X (xCC)                                                 *
 *    - Centroid.Y (yCC)                                                 *
 * Returned data: -the values of the integrals up to 3rd order           *
 *                 i.e. OrderX + OrderY <= 3                             *
 * Obs. Use P1&P2 of the geom. coeff. as OrderX and OrderY               *
 ************************************************************************/
double PolynomLineIntegration(const double & N1x, const double & N1y,
			      const double & N2x, const double & N2y,
			      const double & xCC, const double & yCC,
			      const int OrderX,   const int OrderY){
  
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
