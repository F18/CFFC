/* BGK1DVector.cc:  Source file for BGK1D_Vector class. */


#include"BGK1DVector.h"


/******************************************************//**
 * Static Variables
 **********************************************************/
int BGK1D_Vector::m_length = -1; //default value of -1
                                 //will prohibit creation
                                 //of vectors prior to setting
                                 //the length (which must be 
                                 //constant throughout the simulation).

int BGK1D_Vector::m_set = 0;

double BGK1D_Vector::m_delta_v          = ZERO;
ColumnVector BGK1D_Vector::m_velocities = ColumnVector();

/********************************************************
 * Function: BGK1D_Vector::moment(i)                    *
 *                                                      *
 * return ith velocity moment of distribution           *
 *                                                      *
 ********************************************************/
double BGK1D_Vector::moment(int n) const {
  double moment(0.0);

  for(int i=0; i<m_length; ++i) {
    moment += (*this)[i] * m_delta_v * pow(velocity(i),n);
  }

  return moment;
}

/********************************************************
 * Function: BGK1D_Vector::Maxwell_Bolzmann(rho,u,p)    *
 *                                                      *
 * Creates discretization of MB distribution function   *
 * with density rho, velocity u, and pressure p.        *
 * Note: The density, velocity and pressure of the      *
 *       discrete DF will not be equal to the given     *
 *       parameters.  This is a distribution function   *
 *       which simply contains samples of the           *
 *       continuous one.                                *
 *                                                      *
 ********************************************************/
void BGK1D_Vector::Maxwell_Boltzmann(double rho, double u, double p) {
  double v;
  double B(rho/(2.0*p));

  double AA(-B*u*u+log(rho*sqrt(B/PI)));
  double BB(2.0*B*u);
  double CC(-B);

    for(int i=0; i<m_length; ++i) {
      v = m_velocities[i];
      (*this)[i] = exp (AA + BB*v + CC*v*v);
    }

  return;
}
