#ifndef _LESPREMIXED2D_TURB_INCLUDED
#define _LESPREMIXED2D_TURB_INCLUDED



// Required C++ libraries
#include <cmath> 
#include <complex>
#include <cstdlib>     // defines the drand48() function
#include <ctime>       // defines the time() function
#include <iostream>
#include <fstream>

/*===============================================================*\
   Required header file for Fourier transforms, which are
   performed using the FFTW library. FFTW stands for the Fastest 
   Fourier Transform in the West.
   http://www.fftw.org/
\*===============================================================*/
#ifndef _FFTW_INCLUDED
#include "fftw3.h"
#endif //_FFTW_INCLUDED

#ifndef _MATH_MACROS_INCLUDED
#include "../Math/Math.h"
#endif // _MATH_MACROS_INCLUDED

#ifndef _CFD_INCLUDED
#include "../CFD/CFD.h"
#endif // _CFD_INCLUDED

using namespace std;


// Complex type definition 
typedef  complex<double>  Complex;

// Constants
const Complex  I(0.0, 1.0);      // sqrt(-1.0)




/*======================================*\
            Inline functions                 
\*======================================*/


// x-direction wave number
inline double k_1(const int &n1) {
  return double(n1);
}

inline double k_1(const int &n1, const double &L1) {
  return TWO*PI*double(n1)/L1;
}


// y-direction wave number
inline double k_2(const int &n2) {
  return double(n2);
}

inline double k_2(const int &n2, const double &L2) {
  return TWO*PI*double(n2)/L2;
}


// Random number generator
inline double random_double(){
  return  drand48();
}


// Dissipation of turbulence kinetic energy using a first-order backward
// finite difference
inline double Dissipation_of_TKE(const double &dt, 
				 const double &TKEn, 
				 const double &TKEo) {
  return ((TKEn-TKEo)/dt);
}


// Integral length scale based on the rate of dissipation of turbulence 
//  kinetic energy
inline double Integral_Length_Scale(const double &dissipation, 
				    const double &u_prime) {
  return fabs(u_prime*u_prime*u_prime/dissipation);
}

// Kolmogorov scale
inline double Kolmogorov_Scale(const double &viscosity,
                               const double &dissipation) {
  double vis3;
  vis3 = viscosity*viscosity*viscosity;
  return pow(vis3/dissipation, 0.25);
}


// Prescribed energy spectrum
inline double Energy_Spectrum(const double &abs_wave_num, const int &spectrum_flag){
  double k, kp, kd, eps, u, Lp, EE;
  double C, A, alpha, a_s;
  int s;
 
  k = abs_wave_num;
  
  switch (spectrum_flag) {
    case LEE_REYNOLDS :
    /*****  Lee and Reynolds  *****/
    kp = 4.0, 
      C = 1.0;
    if (k<=kp) {
      EE = C*k*k;
    } else {
      EE = C*kp*kp*pow(k/kp, -5.0/3.0);
    }
    break;

    /*****  Laval and Nazarenko paper  *****/
    case LAVAL_NAZARENKO :
    //double EE, kp, C;
    C = 1.0;
    kp = 4.0;
    EE = C*pow(k, 1.0)*exp(-pow(k/kp, 2.0));
    break;

    /*****   von Karman-Pao   *****/
    case VON_KARMAN_PAO :
    //double A, u, eps, kp, kd, alpha, EE, Lp;
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp = 3.0;   kd = 1000.0 /*362.0*/;
    u = 0.107;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;

    /*****  Haworth and Poinsot paper  *****/
    case HAWORTH_POINSOT :
    //double EE, kp, 
    u = 2.5;  Lp = TWO*PI/6;  // kp=4.0, u=2.5  kp=8
    kp = TWO*PI/Lp;
    EE = (32.0/3.0) * sqrt(2.0/PI)* (u*u/kp) * pow(k/kp, 4.0) * exp(-2.0*(k/kp)*(k/kp));
    break;

    /*****  Chasnov 1996  *****/
    case CHASNOV :  
    //double a_s, EE, kp = 4.0 /*8.0 20.0 4.0*/, u = 0.095; /*0.1 28.3 21.21  0.001*/  
    //int s = 3;
    s = 3;
    kp = 4.0;
    u = 0.009; 
    // u = 0.001  ->  Re_lambda = 10 
    // u = 0.009  ->  Re_lambda = 98
    // u = 0.095  ->  Re_lambda = 950
   
    a_s = pow(2.0*double(s)+1.0, s+1)/(factorial(s)*pow(2.0, s));
    EE = (HALF*a_s*u*u/kp)*pow(k/kp, 2*s+1);
    EE = EE*exp(-(double(s)+HALF)*pow(k/kp, 2));
    break;


    /*****   Bell & Day report   *****/
    case BELL_DAY :
    //  kd = 1/(2*dx)
    kp = 3.0;   kd = ONE/0.576E-3;
    EE = pow(k/kp, 4.0) * exp(-9.0*pow(k/kd, 4.0/3.0)/4.0);
    EE /= pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;


    /*****   von Karman-Pao   *****/
    default :
    //double A, u, eps, kp, kd, alpha, EE, Lp;
    A = 1.5;    alpha = 1.5;
    //Lp = TWO*PI/ONE;
    kp = 3.0;   kd = 1000.0 /*362.0*/;
    u = 0.107;    eps = 2.73E-3; //2.73E-03;
    EE = (A/eps)*pow(u, 5.0);
    EE = EE*pow(k/kp, 4.0) * exp(-3.0*alpha*pow(k/kd, 4.0/3.0)/2.0);
    EE = EE/pow(1.0+(k/kp)*(k/kp), 17.0/6.0);
    break;    

  }  // end switch
  
  if (k==0.0) {
    return 0.0;
  } else {
    return EE;
  }
 
}


// Function taken from Rogallo's report, 1981
inline Complex alpha_Rogallo(const double &abs_wave_num, const double &theta, const int &spectrum_flag){
  double E, k;
  k = abs_wave_num;
  E = Energy_Spectrum(k, spectrum_flag);
 
  if (k == 0) {
    return (0.0, 0.0);
  } else {
    return sqrt(E/(PI*k*k))*exp(I*theta);
    //return sqrt(E/(PI*k))*exp(I*theta);
  }

}


#endif // _LESPREMIXED2D_TURB_INCLUDED
