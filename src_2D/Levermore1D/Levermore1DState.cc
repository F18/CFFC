/* Levermore1DState.cc:  Source file for Levermore1D_XState classes. */

#ifndef _LEVERMORE1D_STATE_INCLUDED
#include "Levermore1DState.h"
#endif //_LEVERMORE1D_STATE_INCLUDED

/********************************************************
 * Static member variables                              *
 ********************************************************/
double Levermore1D_weights::particle_mass = (MOLE_WT_ZB/(AVOGADRO*THOUSAND));

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
 * Function: Levermore1D_cState::find_real_L            *
 *                                                      *
 * For high-order closure, the highest order            *
 * coefficients may be zero when the distribution       *
 * distribution function is actually of lower-order.    *
 * This function finds the "true" length of the weights *
 * vector.                                              *
 *                                                      *
 ********************************************************/
int Levermore1D_cState::find_real_L(const Levermore1D_weights &A) const {
  double moment_test(0.0), min_diff(MILLION*m_values[length-1]);
  int real_L;

  for(int i=3; i <= length; i = i+2) {
    if(A[i] != 0.0) {
      moment_test = moment_series(length-1,A,i);
      if( fabs(m_values[length-1]-moment_test) < min_diff) {
	min_diff = fabs(m_values[length-1]-moment_test);
	real_L = i;
      }
    }
  }
  return real_L;
}

/********************************************************
 * Function: Levermore1D_cState::moment                 *
 *                                                      *
 * Return value of velocity moment              .       *
 *                                                      *
 ********************************************************/
double Levermore1D_cState::moment(int n, const Levermore1D_weights &A) const {
  if(n<length) return m_values[n];
  //else determine real_L (ie if the higher order coefficients are zero or not)
  return moment_series(n,A,find_real_L(A));
}

double Levermore1D_cState::moment(int n, const Levermore1D_weights &A, int real_L) const {
  if(n<length) return m_values[n];
  //else
  return moment_series(n,A,real_L);
}

double Levermore1D_cState::moment_series(int n, const Levermore1D_weights &A, int real_L) const {
  double _moment(0.0); //underscore to differentiate from function name
  double term, max_term(0.0);

  term = (double)(n-real_L+2)*moment(n-real_L+1,A);
  max_term = fabs(term);

  _moment += term;

  for(int i=1; i<real_L-1; ++i) {
    term = (double)(i) * A[i+1]*moment(n-real_L+i+1,A);
    if(fabs(term) > max_term) max_term = fabs(term);
    _moment += term;
  }

  _moment /= (-(double)(real_L-1)*A[real_L]);

  return _moment;
}

/********************************************************
 * Function: Levermore1D_cState::in_sync_with           *
 *                                                      *
 * Use Known relations to determine if U and A are      *
 * close to being syncronized.                          *
 *                                                      *
 ********************************************************/
int Levermore1D_cState::in_sync_with(const Levermore1D_weights &A) const {
  double test1(0.0), test2( (*this)[1] );

  for(int i=1; i < length; ++i) {
    test1 += (double)i * A[i+1] * (*this)[i];
    test2 += (double)i * A[i+1] * (*this)[i+1];
  }

  //scale 'em

  test1 = fabs(test1 * ((*this)[2]/(*this)[1]) );
  test2 = fabs(test2 / (*this)[1] );

  return ( test1 < 1.0e-3 &&
	   test2 < 1.0e-3 && A[length] <= 0.0);
}

/********************************************************
 * Function: Levermore1D_cState::relative_error         *
 *                                                      *
 * Find the relative error between two conserved        *
 * states.                                              *
 *                                                      *
 ********************************************************/
double Levermore1D_cState::relative_error(const Levermore1D_cState &U2) const {
  double error(sqr( ((*this)[1] - U2[1])/(*this)[1] ));
  for(int i = 3; i <= length; i=i+2) {
    error += sqr( ((*this)[i] - U2[i])/(*this)[i] );
    error += sqr( ((*this)[i-1] - U2[i-1])/
		  (*this)[i]*(*this)[1]/(*this)[3]);
  }
  return sqrt(error)/(double)length;
}

/********************************************************
 * Function: Levermore1D_weights::                      *
 *                integrate_conserved_moment            *
 *                                                      *
 * Integrate ith conserved moment                       *
 *                                                      *
 ********************************************************/
double Levermore1D_weights::integrate_conserved_moment(int i) const {
  return integrate_conserved_moment_pos(i) + integrate_conserved_moment_neg(i);
}

double Levermore1D_weights::integrate_conserved_moment_neg(int i) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
  v = -dx/2.0;
  do{
    v += dx;
    temp = velocity_weighted_value_at(v, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  return sum;
}

double Levermore1D_weights::integrate_conserved_moment_pos(int i) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
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
  return integrate_random_moment_pos(i,u) + integrate_random_moment_neg(i,u);
}

double Levermore1D_weights::integrate_random_moment_neg(int i, double u) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
  v = -dx/2.0;
  do{
    v += dx;
    temp = random_velocity_weighted_value_at(v, u, i)*dx;
    sum += temp;
  } while(value_at(v) > 1e-14*dx);

  return sum;
}

double Levermore1D_weights::integrate_random_moment_pos(int i, double u) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = 0.1;
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
int Levermore1D_weights::set_from_U(const Levermore1D_cState &U) {
  if(m_values[length-1] >= 0) {MaxBoltz(U);}

  double rel_err, precon;
  Levermore1D_cState U_temp(*this);
  ColumnVector A_step(Levermore1D_Vector::get_length());
  ColumnVector rhs(Levermore1D_Vector::get_length());
  ColumnVector Plambda(Levermore1D_Vector::get_length());
  DenseMatrix d2hda2(Levermore1D_Vector::get_length(),
			 Levermore1D_Vector::get_length());
  int count(0);

  for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
    rhs[i] = U[i+1]-U_temp[i+1];
  }
  rel_err = U_temp.relative_error(U);

  while( rel_err > 1e-10) { //fix tolerance later
    d2hda2 = U_temp.d2hda2(*this);
//    //precondition d2hda2 & rhs
//    for(int i=0;i<length;++i) {
//      precon = d2hda2( ((i+1)/2)*2 ,0);
//      precon = 1.0;
//      for(int j=0;j<length;++j) {
//    	d2hda2(i,j) /= precon;
//      }
//      rhs(i) /= precon;
//    }
    Solve_LU_Decomposition(d2hda2,rhs,A_step);
    //A_step = rhs;
    //Solve_LS_Householder_F77(d2hda2,A_step,junk,length,length);
    //Solve_LAPACK_dgesv(d2hda2,A_step);
    *this += A_step;
    for(int L = length; L > 1; L = L-2) {
      if(L==3 && (*this)[L] > 0.0 ) return 1;
      if( (*this)[L] < 0.0) break;
      //else
      (*this)[L] = 0.0; (*this)[L-1] = 0.0;
    }
    U_temp.set_from_A(*this);
    rel_err = U_temp.relative_error(U);

    for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
      rhs[i] = U[i+1]-U_temp[i+1];
    }
    ++count;
    if(count > 50) return 1;
  }

  return 0;
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

  int real_L;
  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  real_L = find_real_L(A);

  for(int i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
      dm(i,j) = moment(i+j,A,real_L);
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

  int real_L;
  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  real_L = find_real_L(A);

  for(int i = 0; i < Levermore1D_Vector::get_length(); ++i) {
    for(int j = 0; j < Levermore1D_Vector::get_length(); ++j) {
      dm(i,j) = moment(i+j+1,A,real_L);
    }
  }

  return dm;
}

/********************************************************
 * Function: FluxHLLE                                   *
 *                                                      *
 * Calculate HLLE flux function.                        *
 *                                                      *
 ********************************************************/
Levermore1D_Vector FluxHLLE(const Levermore1D_cState &Ul,
			    const Levermore1D_weights &Al,
			    const double &wavespeed_l,
			    const Levermore1D_cState &Ur,
			    const Levermore1D_weights &Ar,
			    const double &wavespeed_r) {
  Levermore1D_Vector Flux;

  if (wavespeed_l >= ZERO) {
    Flux = Ul.F(Al);
  } else if (wavespeed_r <= ZERO) {
    Flux = Ur.F(Ar);
  } else {
    Flux =   ((Ul.F(Al)*wavespeed_r-Ur.F(Ar)*wavespeed_l)
	      +(Ur-Ul)*(wavespeed_l*wavespeed_r))/
             (wavespeed_r-wavespeed_l);
  } /* endif */

  return Flux;
}

/********************************************************
 * Function: FluxKinetic                                *
 *                                                      *
 * Calculate "Kinetic" flux function.                   *
 *                                                      *
 ********************************************************/
Levermore1D_Vector FluxKinetic(const Levermore1D_weights &Al,
			       const Levermore1D_weights &Ar) {
  Levermore1D_Vector Flux;

  for(int i=1; i<=Levermore1D_Vector::get_length(); ++i) {
    Flux[i] = Al.integrate_conserved_moment_pos(i);
    Flux[i] = Ar.integrate_conserved_moment_neg(i);
  }

  return Flux;
}
