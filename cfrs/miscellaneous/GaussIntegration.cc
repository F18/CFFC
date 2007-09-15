//////////////////////////////////////////////////////////////////////////
// This program codes an example of using the Gauss-Legendre
// method to estimate the integral 
//
//          I = int[f(x),x=a..b],
//
// where f(x) is a user-defined function.  For this example
// we will assume f(x) is e^x and a=0 b=1.
//
//                     thomas l. gibson  April 24, 2000
//
//          updated by thomas l. gibson  April 20, 2004
//////////////////////////////////////////////////////////////////////////

#include <cstdlib>               // Needed for system calls
#include <iostream>              // Needed for cout and cin
#include <fstream>               // Needed for file output
#include <iomanip>               // Needed for setw function in file out
#include <cmath>                 // Needed for math functions
#include <cassert>               // Needed for assert
using namespace std;

double f(double t);              // Function prototype for integrand
// G-L pts & wts function 
int gauss_legendre(double , double , double [], double [], int); 

const int nmax = 128;            // Dimension for 128-pt G-L quadrature

int main()
{
  double a,b,error,relerr,sum,x,exact;
  double xx[nmax],ww[nmax];      // Arrays for G-L pts & wts
  int n,i,doit,ierr;
  char ask[2];
  exact=exp(1.0) - 1.0;          // Exact value of integral
  doit=1;
//   a=0.0;                         // Lower limit of integration
//   b=1.0;                         // Upper limit of integration
  a = -125.0;
  b = 1250.0;
  fstream prt_1;
  prt_1.open("gauss_results.dat", ios::out);
  prt_1 << "Integration results using Gauss-Legendre quadrature \n";
  prt_1 << "\t for the integral int[e^x,x=0..1] = "
        << exact << "\n";
  system("clear");
  while(doit)
  {
    cout << "Input the number of G-L quadrature points (0 < n < "
         << nmax + 1 << "): \n";
    cout << "\t n = ? ";
    cin >> n;
    if((n <= 0) || (n > nmax))
    {
      cout << "Bad News! Your choice is illegal.\n";
      exit(0);
    }
    cout << "\n" << "We are using " << n <<" G-L points \n";

    ierr = gauss_legendre(-1.0,1.0,xx,ww,n); // Compute G-L pts & wts
    sum=0.0;
    for (i=0; i<n; i++)
    {
      x = 0.5*(b + a + (b - a)*xx[i]);      // transform for pts
      sum += ww[i] * f(x);
    }
    sum *= 0.5*(b - a);                     // transform for wts
    error=sum - exact;
    relerr=(error/exact)*100.0;
    cout << "For " << n << " G-L points, "
         << " Int = " << sum << " error = " << error
         << " %rel. error = " << relerr << "\n";
    prt_1 << "n = " << setw(5) << n 
          << " Int = " << setw(10) << sum 
          << " error = " << setw(12) << error
          << " %rel. error = " << setw(12) << relerr << "\n";
    cout << "\n \n" << "Change the number of points and try again? (y/n) ";
    cin >> ask;
    if(ask[0] == 'n' || ask[0] == 'N')
    {
    doit=0;
    }
  }
  return 0;
}
//////////////////////////////////////////////////////////////////////////   
// This function calculates the integrand
// for f(x)=e^x
//               thomas l. gibson  April 24, 2000
//////////////////////////////////////////////////////////////////////////
double f(double t)
{
  double y;
  //  y=exp(t);
  y = pow(t,6)-4*pow(t,4)+2*t;
  return y;
}
//////////////////////////////////////////////////////////////////////////
//  This function computes the points and weights for 
//  an n-point gauss-legendre quadrature.  This is based
//  on the routine gauleg in "Numerical Recipes in C," by
//  W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
//  B.P. Flannery.
//  
//  This function needs <math.h>
//  N.B. The points x[] and weights w[] arrays begin at 
//  i = 0.
//
//              thomas l. gibson  April 22, 2000
//////////////////////////////////////////////////////////////////////////

int gauss_legendre(double x1, double x2, double x[], double w[], int n)
{
  double z1,z,xm,xl,pp,p3,p2,p1;
  static const double pi = M_PI;
  static const double EPS = 1.0e-16;    // relative precision for 
                                        // the Newton-Raphson loop
  int m;
  int ierr = -1;                        // return status
  register int i,j;

  m = (n + 1)/2;                        // roots are symmetric - find half
  xm = 0.5*(x2 + x1);
  xl = 0.5*(x2 - x1);
  for (i=1; i <= m; i++)
  {
    int i_guard = 1;                    // added by tlg  April 22, 2004
    z = cos(pi * (i - 0.25)/(n + 0.5)); // approximation to ith root
    do                                  // use Newton-Raphson to find roots
    {
      p1 = 1.0;
      p2 = 0.0;
      for (j = 1; j <= n; j++)         
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3)/j; // recurrence relation
      }
      pp = n * (z * p1 - p2)/(z * z - 1.0);           // derivative of p_n
      z1 = z;
      z = z1 - p1/pp;                  // Newton-Raphson step
      i_guard++;                       // Guard value to end loop added
                                       // by tlg April 22, 2004
    } while ( (fabs(z - z1) > EPS) && (i_guard < 100) );//tlg April 22, 2004
    x[i - 1] = xm - xl * z;        // pts & wts arrays
    cout << "x[" << i-1 << "]=" << setprecision(20) << x[i-1] << "\t";
    x[n - i] = xm + xl * z;        // start at i = 0
    w[i - 1] = 2.0 * xl/((1.0 - z*z) * pp * pp);
    cout << "w[" << i-1 << "]=" << w[i-1] << endl;
    cout << "x[" << n-i << "]=" << x[n-i] << endl;
    w[n - i] = w[i - 1];
  }
  ierr = 0;
  return ierr;
}
//////////////////////////////////////////////////////////////////////////
