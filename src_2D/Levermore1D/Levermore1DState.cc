/* Levermore1DState.cc:  Source file for Levermore1D_XState classes. */

#ifndef _LEVERMORE1D_STATE_INCLUDED
#include "Levermore1DState.h"
#endif //_LEVERMORE1D_STATE_INCLUDED

/********************************************************
 * Static member variables                              *
 ********************************************************/
double Levermore1D_weights::particle_mass = 1.0e-20; //something small?

/********************************************************
 * Function: Levermore1D_pState::set_from_U             *
 *                                                      *
 * Convert a converved vector into a primitive one.     *
 *                                                      *
 ********************************************************/
void Levermore1D_pState::set_from_U(const Levermore1D_cState &U) {

  for(int i=0; i<length; ++i) {
    if(i>1) {
      m_values[i] = U[i+1]-conserved_extras(i);
    } else if(i==1) {
      m_values[i] = U[2]/U[1];
    } else {  //i=0
      m_values[i] = U[1];
    }
  }

  return;
}

/********************************************************
 * Function: Levermore1D_cState::set_from_W             *
 *                                                      *
 * Convert a primitive vector into a conserved one.     *
 *                                                      *
 ********************************************************/
void Levermore1D_cState::set_from_W(const Levermore1D_pState &W) {

  for(int i=0; i<length; ++i) {
    if(i>1) {
      m_values[i] = W[i+1]+W.conserved_extras(i);
    } else if(i==1) {
      m_values[i] = W[1]*W[2];
    } else {  //i=0
      m_values[i] = W[1];
    }
  }

  return;
}

/********************************************************
 * Function: Levermore1D_pState::conserved_extras       *
 *                                                      *
 * Find the difference between the "i" entry in the     *
 * primitive and conserved vectors.  This requires all  *
 * previous entries in the primitive vector be up to    *
 * date.                                                *
 *                                                      *
 ********************************************************/
double Levermore1D_pState::conserved_extras(int i) const {
  int pf(1);
  return m_values[1]*conserved_extras_recursive(i-1,pf,i,1);
}

/********************************************************
 * Function: Levermore1D_pState::                       *
 *              conserved_extras_recursive              *
 *                                                      *
 * Used for the conserved_extras function.              *
 *                                                      *
 ********************************************************/
double Levermore1D_pState::conserved_extras_recursive(int i, int &pf, int pf_num, int pf_den) const {
  if(i==1) {
    return m_values[1]*m_values[0];
  } else {
    pf = (pf*pf_num)/pf_den;
    return (double)pf * m_values[i]  + m_values[1] * conserved_extras_recursive(i-1,pf,pf_num-1,pf_den+1);
  }
}

/********************************************************
 * Function: Levermore1D_cState::set_from_A             *
 *                                                      *
 * Convert a weigths vector into a conserved one.       *
 *                                                      *
 ********************************************************/
void Levermore1D_cState::set_from_A(const Levermore1D_weights &A) {
  for(int i=0; i<length; ++i) {
    m_values[i] = A.integrate_conserved_moment(i);
  }
  return;
}

/********************************************************
 * Function: Levermore1D_pState::set_from_A             *
 *                                                      *
 * Convert a weigths vector into a primitive one.       *
 *                                                      *
 ********************************************************/
void Levermore1D_pState::set_from_A(const Levermore1D_weights &A) {


  m_values[0] = A.integrate_conserved_moment(0);
  m_values[1] = A.integrate_conserved_moment(1)/m_values[0];

  for(int i=2; i<length; ++i) {
    m_values[i] = A.integrate_random_moment(i,m_values[1]);
  }
  return;
}

/********************************************************
 * Function: Levermore1D_cState::moment                 *
 *                                                      *
 * Return value of velocity moment              .       *
 *                                                      *
 ********************************************************/
double Levermore1D_cState::moment(int n, const Levermore1D_weights &A) const {
  if(n<length) return m_values[n];

  double _moment(0.0); //underscore to differentiate from function name
  int L(length);

  //while(fabs(A[L])<1e-6) {L -=2;}  //There must be something better than this.

  _moment += (double)(n-L+2)*moment(n-L+1,A);

  for(int i=1; i<L-1; ++i) {
    _moment += (double)(i) * A[i+1]*moment(n-L+i+1,A);
  }

  _moment /= -(double)(L-1)*A[L];
  return _moment;
}

/********************************************************
 * Function: Levermore1D_weights::                      *
 *                integrate_conserved_moment            *
 *                                                      *
 * Integrate ith conserved moment                       *
 *                                                      *
 ********************************************************/
double Levermore1D_weights::integrate_conserved_moment(int i) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
  v = -dx/2.0;

  do{
    v += dx;
    temp = velocity_weighted_value_at(v, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  v = dx/2.0;
  do{
    v -= dx;
    temp = velocity_weighted_value_at(v, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  return sum;
}

/********************************************************
 * Function: Levermore1D_weights::                      *
 *                integrate_random_moment               *
 *                                                      *
 * Integrate ith random moment                          *
 *                                                      *
 ********************************************************/
double Levermore1D_weights::integrate_random_moment(int i, double u) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
  v = -dx/2.0;

  do{
    v += dx;
    temp = random_velocity_weighted_value_at(v, u, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  v = dx/2.0;
  do{
    v -= dx;
    temp = random_velocity_weighted_value_at(v, u, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  return sum;
}


/********************************************************
 * Function: Levermore1D_weights::set_from_U            *
 *                                                      *
 * Convert a conserved vector into weights.             *
 *                                                      *
 ********************************************************/
void Levermore1D_weights::set_from_U(const Levermore1D_cState &U) {

  Levermore1D_cState U_temp(*this);
  Levermore1D_weights A_step;
  ColumnVector rhs(Levermore1D_Vector::get_length());
  DenseMatrix d2hda2_inv(Levermore1D_Vector::get_length(),
			 Levermore1D_Vector::get_length());

  for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
    rhs[i] = U[i+1]-U_temp[i+1];
  }

  while(fabs(rhs[1]) > 1e-10) { //fix tolerance later
    d2hda2_inv =U_temp.d2hda2(*this).pseudo_inverse();
    A_step = d2hda2_inv*rhs;
    *this += A_step;
    while( (*this).m_values[get_length()-1] > 0 ) {
      A_step = A_step*0.5;
      *this -= A_step;
    }
    U_temp.set_from_A(*this);
    for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
      rhs[i] = U[i+1]-U_temp[i+1];
    }
  }

  return;
}

/********************************************************
 * Function: Levermore1D_weights::set_from_W            *
 *                                                      *
 * Convert a primitive vector into weights.             *
 *                                                      *
 ********************************************************/
void Levermore1D_weights::set_from_W(const Levermore1D_pState &W) {
  Levermore1D_cState U(W);
  set_from_U(U);
  return;
}

/********************************************************
 * Function: Levermore1D_cState::d2hda2                 *
 *                                                      *
 * Calculate Hessian of density potential.              *
 *                                                      *
 ********************************************************/
DenseMatrix Levermore1D_cState::d2hda2(const Levermore1D_weights &A) const {

  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  for(int i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
      dm(i,j) = moment(i+j,A);
    }
  }
  return dm;
}

/********************************************************
 * Function: Levermore1D_cState::d2jda2                 *
 *                                                      *
 * Calculate Hessian of flux potential.                 *
 *                                                      *
 ********************************************************/
DenseMatrix Levermore1D_cState::d2jda2(const Levermore1D_weights &A) const {

  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  for(int i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
      dm(i,j) = moment(i+j+1,A);
    }
  }
  return dm;
}

