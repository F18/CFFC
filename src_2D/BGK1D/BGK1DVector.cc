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
  ColumnVector V(MB_coefs(rho,u,p));
  fill_with_MB(V(0), V(1), V(2));
  return;
}

ColumnVector BGK1D_Vector::MB_coefs(double rho, double u, double p) const {
  double B(rho/(2.0*p));
  ColumnVector V(3);

  V(0) = -B*u*u+log(rho*sqrt(B/PI));
  V(1) = 2.0*B*u;
  V(2) = -B;

  return V;
}

void BGK1D_Vector::fill_with_MB(double AA, double BB, double CC) {
  double v;
  for(int i=0; i<m_length; ++i) {
    v = m_velocities[i];
    (*this)[i] = exp (AA + BB*v + CC*v*v);
  }
  return;
}

/********************************************************
 * Function: BGK1D_Vector::discrete_Maxwell_Bolzmann    *
 *                                                      *
 * Creates a discretization which has density rho,      *
 * velocity u, and pressure p and maximizes entropy.    *
 * Note: It is a Gaussian distribution function but     *
 *       has different coefficients than the continuous *
 *       case.                                          *
 *                                                      *
 ********************************************************/
int BGK1D_Vector::discrete_Maxwell_Boltzmann(const BGK1D_Vector &V_in) {
  double rho, u, p;
  rho = V_in.moment(0);
  u = V_in.moment(1)/rho;
  p = V_in.moment(2) - rho * u * u;
  return discrete_Maxwell_Boltzmann(rho, u, p);
}

int BGK1D_Vector::discrete_Maxwell_Boltzmann(double rho, double u, double p) {

  int count(0), j_max, j_min;
  double term;
  ColumnVector coefs(3), want(3), have(3), rhs(3), update(3);
  DenseMatrix  lhs(3,3);

  coefs = MB_coefs(rho,u,p);
  fill_with_MB(coefs(0), coefs(1), coefs(2));

  want[0] = rho;
  want[1] = rho*u;
  want[2] = p + rho*u*u;

  have[0] = moment(0);
  have[1] = moment(1);
  have[2] = moment(2);

  while( fabs( (want[0]-have[0])/want[0] ) > tolerance() ||
	 fabs( (want[1]-have[1])/want[1] ) > tolerance() ||
	 fabs( (want[2]-have[2])/want[2] ) > tolerance() ) {

    rhs = want-have;

    for(int i = 0; i <= 4; ++i) {
      term = moment(i);
      j_min = max(0,i-2);
      j_max = min(i,2);
      for(int j = j_min; j <= j_max; ++j) lhs(i-j,j) = term;
    }

    Solve_LU_Decomposition(lhs,rhs,update);

    coefs += update;
    fill_with_MB(coefs(0), coefs(1), coefs(2));

    have[0] = moment(0);
    have[1] = moment(1);
    have[2] = moment(2);

    ++count;
    if(count > 100) {
      cout << endl << endl << "Error, could not create discrete distribution function" << endl
	   << "with a pressure of " << p << " Pa, a bulk velocity of" << endl
	   << u << " m/s, and a density of " << rho << " Kg/m^3." << endl << endl;
      return 1;
    }

  }

  return 0;
}
