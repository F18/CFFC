/******************************************************
 ******************************************************
 ** Program:  simpson.h                              **
 **                                                  **
 ** Description:  This program uses adaptive         **
 **               simpson quadrature to integrate    **
 **               functions numerically.  See:       **
 **                                                  **
 ** W. Gander and W. Gautschi, "Adaptive Quadrature -**
 ** Revisisted." BIT, Vol. 40, 2000, pp. 84-101.     **
 ** This document is also available at               **
 ** http://www.inf.ethz.ch/personal/gander           **
 **                                                  **
 ** Author: Marc "T-Bone" Charest                    **
 **                                                  **
 ** Revision:  Date       Initials  Change           **
 **            23/04/2007 MRC       Original creation**
 **                                                  **
 ******************************************************
 ******************************************************/

#ifndef _SIMPSON_INCLUDED
#define _SIMPSON_INCLUDED


//--------------------------------------------------//
// SAMPLE USAGE
//--------------------------------------------------//
// double fnc(int ndim, double *x, void *params)
// {
//   if (x[0]<1) return x[0]+1.0;
//   else if (x[0]>=1 && x[0]<=3) return 3.0-x[0];
//   else return 2.0;
  
// }

// int main() {
  
//   // declares
//   simp_function F;
//   simp_state S;
//   simp_params P;
//   int fevals;
//   double val;
//   int ret;

//   // alocate memory
//   malloc_simp_struc( 4, F, S );
//   init_simp_struc( F, S );
  
//   // setup function
//   F.f = fnc;
//   F.xmin[0] = 0; F.xmax[0] = 5;
//   F.xmin[1] = 0; F.xmax[1] = 1;
//   F.xmin[2] = 0; F.xmax[2] = 1;
//   F.xmin[3] = 0; F.xmax[3] = 1;

//   // setup integration parameters
//   P.maxevals = 100000;
//   P.tol = 1e-6;

//   // perform integral
//   ret = adaptsim( F, S, P, fevals, val );

//   // output
//   cout.setf(ios::fixed);
//   cout.precision(14);
//   cout.setf(ios::left);
//   cout << "Result          : " << setw(22) << val << endl; 
//   cout << "Function Evals. : " << fevals << endl;
//   cout << "Flag            : " << ret << endl;

//   // free memory
//   free_simp_struc( F, S );

//   return 0;

// }
//--------------------------------------------------//
// END SAMPLE USAGE
//--------------------------------------------------//





#include <math.h>
#include <iostream>
#include <iomanip>
using namespace std;


/* Constants */

const int SIMPSON_SUCCESS      = 0;  // converged successfully
const int SIMPSON_MACHINE_ZERO = 1;  // cannot refine any further
const int SIMPSON_MAX_ITER     = 2;  // max iteration count exceeded
const int SIMPSON_NAN          = 3;  // NaN



/* typedefs */

typedef double (*integrand) (int ndim, double *x, void *params);


/* structs */

struct simp_function {
  integrand f;
  int ndim;
  double *xmin;
  double *xmax;
  void *params;
};

struct simp_state {
  int dim;
  double *x;
};

struct simp_params {
  int maxevals;
  double tol;
};




/* Function Prototypes */

// machine epsilon
double epsilon( const double d ); 


// main N-Dim outer control loop for simpson integration
int adaptsim( const simp_function &F,
	      simp_state &S,
	      const simp_params &P,
	      int &feval,
	      double &val );

// main workhorse for 1D simpson integration
int adaptsimstp( const simp_function &F,
		 simp_state &S,
		 const double a, const double b, 
		 const double fa, const double fm, 
		 const double fb, const double is,
		 const double hmin, 
		 const simp_params &P,
		 int &fevals,
		 double &val );

// recursive function evaluation
int evaluate( const simp_function &F,
	      simp_state &S,
	      const simp_params &P,
	      int &fevals,
	      double &val );

// allocate and deallocate structs
void malloc_simp_struc( const int dim, simp_function &F, simp_state &S );
void free_simp_struc( simp_function &F, simp_state &S );
void init_simp_struc( simp_function &F, simp_state &S );
   


/****************************************************
 * Allocate and free memory for necessary structs   *
 ****************************************************/
inline void malloc_simp_struc( const int dim, 
			       simp_function &F, 
			       simp_state &S ) {
  F.ndim = dim;
  F.xmin = new double[F.ndim]; 
  F.xmax = new double[F.ndim]; 
  S.x    = new double[F.ndim]; 
}


inline void free_simp_struc( simp_function &F, simp_state &S ) 
{ delete[] F.xmin;  delete[] F.xmax;  delete[] S.x; }


inline void init_simp_struc( simp_function &F, simp_state &S )
{
  S.dim = 0;
  for ( int i=0; i<F.ndim; i++) {
    F.xmin[i] = 0;
    F.xmax[i] = 0;
    S.x[i]    = 0;
  }
}


/****************************************************
 * Adaptive Simpson integration function            *
 * (main workhorse)                                 *
 ****************************************************/
inline int adaptsimstp( const simp_function &F,
			simp_state &S,
			const double a, const double b, 
			const double fa, const double fm, 
			const double fb, const double is,
			const double hmin, 
			const simp_params &P,
			int &fevals,
			double &val )
{
  
  //declares
  double m, h, fml, fmr;
  double i1, i2;
  int which_dim = S.dim;
  int err;
  double temp;
  
  // compute the midpoint and the panel size
  m = (a+b)/2.0;
  h = (b-a)/4.0;
  
  // evaluate the function at the panel midpoints
  S.x[which_dim] = a+h;  S.dim = which_dim;  
  err = evaluate( F, S, P, fevals, fml );
  S.x[which_dim] = b-h;  S.dim = which_dim;  
  err = evaluate( F, S, P, fevals, fmr );
  fevals += 2;
  
  // use 3pt Simsons  rule
  i1 = (h/1.5) * (fa + 4.0*fm + fb);
  
  // use 5pt Simsons  rule
  i2 = (h/3.0) * (fa + 4.0*fml + 2.0*fm + 4.0*fmr + fb);
  
  // one step of Romberg extrapolation
  i1 = (16.0*i2 - i1)/15.0;
  
  // check NaN
  if ( i1!=i1 || i2!=i2) {
#ifdef _DEBUG_
    cerr << "Error: Infinite or Not-a-Number function value encountered."
	 << endl;
#endif
    val = i1;
    err = SIMPSON_NAN;
  }


  // if convergence is met or we cannot refine any further
  if ( (fabs(i1-i2)<P.tol) || (m<=a) || (b<=m) || (h<hmin) ||
       (fevals>=P.maxevals)  ) {

    // and we cannot refine any further
    if ( (m<=a) || (b<=m) || h<hmin )  {
#ifdef _DEBUG_
      cerr << "Warning: Interval contains no more "
	   << "machine number.  Required tolerance "
	   << "may not be met."
	   << endl;
#endif
      err = SIMPSON_MACHINE_ZERO;
    } 

    // iteration count exceeded
    else if ( fevals>=P.maxevals ) {
#ifdef _DEBUG_
      cerr << "Warning: Maximum number of iterations exceeded."
	   << endl;
#endif
      err = SIMPSON_MAX_ITER;
    } 

    // must be ok
    else {
      err = SIMPSON_SUCCESS;
    }
    /* endif */
    
    
    // set the refined value
    val = i1;
      
    
  // else convergence is still not met, break the interval
  // into 2 and recompute
  } else {
    S.dim = which_dim; 
    err = adaptsimstp( F, S, a, m, fa, fml, fm, is, hmin, P, fevals, temp );
    val = temp;
    S.dim = which_dim; 
    err = adaptsimstp( F, S, m, b, fm, fmr, fb, is, hmin, P, fevals, temp );
    val += temp;
  
  } /* endif */
  
  return err;
  
}



/****************************************************
 * Adaptive simpson integration (control loop)      *
 ****************************************************/
inline int adaptsim( const simp_function &F,
		     simp_state &S,
		     const simp_params &P,
		     int &fevals,
		     double &val )
{

  // declares
  static double monte_carlo[5] = {0.9501, 0.2311, 0.6068, 0.4860, 0.8913};
  double m, fa, fb, fm, ft;
  double y_sum, yy_sum;
  double is;
  int j;
  double hmin;
  int err;
  int which_dim = S.dim;
  double a = F.xmin[which_dim];
  double b = F.xmax[which_dim];
  
  // if this is the first call
  if (which_dim == 0) { val = 0; fevals = 0; }

  // compute the midpoint
  m = (a+b)/2.0;
  
  // compute the function values
  S.x[which_dim] = a;  S.dim = which_dim; 
  err = evaluate( F, S, P, fevals, fa );
  S.x[which_dim] = m;  S.dim = which_dim; 
  err = evaluate( F, S, P, fevals, fm );
  S.x[which_dim] = b;  S.dim = which_dim; 
  err = evaluate( F, S, P, fevals, fb );


  // we need an estimate for the modulus of the integral 
  // here we use a Monte-Carlo estimate which also uses 
  // the function values in the middle and at the end 
  // points of the interval
  y_sum = fa + fm + fb;
  for (yy_sum=0, j=0; j<5; j++) {
    S.dim = which_dim; 
    S.x[which_dim] = a + monte_carlo[j]*(b-a);
    err = evaluate( F, S, P, fevals, ft );
    yy_sum += ft;
  }
  is = (b-a)*(y_sum+yy_sum)/8.0;
  fevals += 8;
  
  // if our estimate is zero, use b-a
  if (is==0) is = b-a;
  
  // determine minimum step size
  hmin = epsilon(b-a);
  
  // call the core integral routine and return the result
  S.dim = which_dim; 
  err = adaptsimstp( F, S, a, b, fa, fm, fb, is, hmin, P, fevals, val );

  return err;

}


/****************************************************
 * Evaluate function                                *
 ****************************************************/
inline int evaluate( const simp_function &F,
		     simp_state &S,
		     const simp_params &P,
		     int &fevals,
		     double &val )
{

  if (S.dim<F.ndim-1) {
    S.dim += 1;
    return adaptsim( F, S, P, fevals, val );
  } else {
    val = F.f(F.ndim, S.x, F.params);
    return SIMPSON_SUCCESS;
  }
}


/****************************************************
 * Determine machine epsilon                        *
 ****************************************************/
inline double epsilon( const double d ) {

  double temp1, temp2, eps;
  temp1 = 1.0;
  
  do {
    eps = temp1;
    temp1 /= 2.0;
    temp2 = d + temp1;
  } while (temp2 > d);

  return eps;
}

#endif // _SIMPSON_INCLUDED
