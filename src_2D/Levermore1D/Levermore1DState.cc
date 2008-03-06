/* Levermore1DState.cc:  Source file for Levermore1D_XState classes. */

#ifndef _LEVERMORE1D_STATE_INCLUDED
#include "Levermore1DState.h"
#endif //_LEVERMORE1D_STATE_INCLUDED

#ifndef _LEVERMORE1D_NUMERICALLIBRARY_WRAPPER_H
#include "Levermore1D_NLwrapper.h"
#endif //_LEVERMORE1D_NUMERICALLIBRARY_WRAPPER_H

#ifndef _NUMERICAL_LIBRARY_INCLUDED
#include "../Math/NumericalLibrary.h"
#endif // _NUMERICAL_LIBRARY_INCLUDED

#define INTEGRATION_RESOLUTION  10.0

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
void Levermore1D_cState::set_from_A(const Levermore1D_weights &A, double us) {
  for(int i=0; i<length; ++i) {
    m_values[i] = A.integrate_conserved_moment(i, us);
  }
  return;
}

/********************************************************
 * Function: Levermore1D_pState::set_from_A             *
 *                                                      *
 * Convert a weigths vector into a primitive one.       *
 *                                                      *
 ********************************************************/
void Levermore1D_pState::set_from_A(const Levermore1D_weights &A, double us) {

  m_values[0] = A.integrate_conserved_moment(0,us);
  m_values[1] = A.integrate_conserved_moment(1,us)/m_values[0];

  for(int i=2; i<length; ++i) {
    m_values[i] = A.integrate_random_moment(i,m_values[1],us);
  }
  return;
}

/********************************************************
 * Function: Levermore1D_cState::moment                 *
 *                                                      *
 * Return value of velocity moment              .       *
 *                                                      *
 ********************************************************/
double Levermore1D_cState::moment(int n, const Levermore1D_weights &A, const double &us) const {
  if(n<length) return m_values[n];
  //return A.integrate_conserved_moment(n, us);
  if(n==length) return moment_series_L(A,us);
  return moment_series(n,A,us);
}

double Levermore1D_cState::moment_series(int n, const Levermore1D_weights &A, const double &us) const {

  double _moment(0.0); //underscore to differentiate from function name
  int i, sign(1);

  _moment = (double)(n-length)*moment(n-length-1,A,us);

  for(i=0; i<length-1; ++i) {
    _moment += ( (double)(i+1)*A[i+2]+
		 sign*STABILIZATION*(double)(length+1)*Pascals_Triangle(length,i)*pow(us,length-i) ) * moment(n-length+i,A,us);
    sign *= -1;
  }
  
  _moment += ( sign*STABILIZATION*(double)(length+1)*Pascals_Triangle(length,i)*pow(us,length-i) ) * moment(n-length+i,A,us);

  _moment /= ((double)(length+1)*STABILIZATION);

  return _moment;
}

double Levermore1D_cState::moment_series_L(const Levermore1D_weights &A, const double &us) const {
  double _moment(0.0);
  int i, sign(1);

  for(i=0; i<length-1; ++i) {
    _moment  += ( (double)(i+1)*A[i+2] +
		  (double)sign*STABILIZATION*(double)(length+1)*Pascals_Triangle(length,i)*pow(us,length-i) ) * moment(i,A,us);
    sign *= -1;
  }

  _moment += ( sign*STABILIZATION*(double)(length+1)*Pascals_Triangle(length,i)*pow(us,length-i) ) * moment(i,A,us);

  _moment /= ((double)(length+1)*STABILIZATION);

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
  double moment1, moment2, test;
  double us(m_values[1]/m_values[0]);
  moment1 = moment_series(length,A,us);
  moment2 = A.integrate_conserved_moment(length,us);
  test = fabs((moment1-moment2)*A[length]);
//  cout << endl << "********************" << endl
//       << (*this) << endl << Levermore1D_cState(A,m_values[1]/m_values[0]) << endl;
//  cout << A << endl;
//  cout << moment1 << " " << moment2 << endl;
//  cout << test << endl; 
//  cout << "Next moment = " << moment_series(length,A,us) << "   " << A.integrate_conserved_moment(length,us) << endl;
//  cout.flush();
  return (test < 1.0e-1);
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
double Levermore1D_weights::integrate_conserved_moment(int i, double us) const {
  return integrate_conserved_moment_pos(i,us) + integrate_conserved_moment_neg(i,us);
}

double Levermore1D_weights::integrate_conserved_moment_pos(int i, double us) const {
//
//  double sum(0.0), term, v_min(0.0), v_max(100.0);
//  Levermore1D_Wrapper Wrapper(this,&Levermore1D_weights::velocity_weighted_value_at,i,us);
//
//  do{
//    term = AdaptiveGaussianQuadrature(Wrapper,v_min,v_max,term,14);
//    sum += term;
//    v_min = v_max;
//    v_max += 100.0;
//  }while(fabs(term)>=1e-14*fabs(sum) && fabs(v_max) < 20.0*100.0);
//
//  if(fabs(v_max) >= 20.0*1000.0) {
//    cout << endl << "Error #121312412" << endl;
//    assert(0==1);
//  }
//
//  return sum;
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = INTEGRATION_RESOLUTION;
  v = -dx/2.0;
  do{
    v += dx;
    temp = velocity_weighted_value_at(v, i, us)*dx;
    sum += temp;
    if(fabs(v)>100000.0) {
      cout << "Error in integration:" << endl
	   << "v = " << v << endl
	   << "value_at = " << value_at(v, us) << endl
	   << "A = " << (*this) << endl;
      assert(1==0);
    }
  } while(value_at(v, us) > 1e-24*dx);

  return sum;
}

double Levermore1D_weights::integrate_conserved_moment_neg(int i, double us) const {
//
//  double sum(0.0), term, v_min(-100.0), v_max(0.0);
//  Levermore1D_Wrapper Wrapper(this,&Levermore1D_weights::velocity_weighted_value_at,i,us);
//
//  do{
//    term = AdaptiveGaussianQuadrature(Wrapper,v_min,v_max,term,14);
//    sum += term;
//    v_max = v_min;
//    v_min -= 100.0;
//  }while(fabs(term)>=1e-14*fabs(sum) && fabs(v_max) < 20.0*100.0);
//
//  if(fabs(v_max) >= 20.0*1000.0) {
//    cout << endl << "Error #717132" << endl;
//    assert(0==1);
//  }
//
//  return sum;
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = INTEGRATION_RESOLUTION;
  v = dx/2.0;
  do{
    v -= dx;
    temp = velocity_weighted_value_at(v, i, us)*dx;
    sum += temp;
    if(fabs(v)>100000.0) {
      cout << "Error in integration:" << endl
	   << "v = " << v << endl
	   << "value_at = " << value_at(v, us) << endl
	   << "A = " << (*this) << endl;
      assert(1==0);
    }
  } while(value_at(v, us) > 1e-24*dx);

  return sum;
}

/********************************************************
 * Function: Levermore1D_weights::                      *
 *                integrate_random_moment               *
 *                                                      *
 * Integrate ith random moment                          *
 *                                                      *
 ********************************************************/
double Levermore1D_weights::integrate_random_moment(int i, double u, double us) const {
  return integrate_random_moment_pos(i, u, us) + integrate_random_moment_neg(i, u, us);
}

double Levermore1D_weights::integrate_random_moment_pos(int i, double u, double us) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = INTEGRATION_RESOLUTION;
  v = -dx/2.0;
  do{
    v += dx;
    temp = random_velocity_weighted_value_at(v, u, i, us)*dx;
    sum += temp;
  } while(value_at(v, us) > 1e-14*dx);

  return sum;
}

double Levermore1D_weights::integrate_random_moment_neg(int i, double u, double us) const {
  double temp(0.0), sum(0.0), pos(0.0), v(0.0);
  double dx = INTEGRATION_RESOLUTION;
  v = dx/2.0;
  do{
    v -= dx;
    temp = random_velocity_weighted_value_at(v, u, i, us)*dx;
    sum += temp;
  } while(value_at(v, us) > 1e-14*dx);

  return sum;
}


/********************************************************
 * Function: Levermore1D_weights::set_from_U            *
 *                                                      *
 * Convert a conserved vector into weights.             *
 *                                                      *
 ********************************************************/
int Levermore1D_weights::set_from_U(const Levermore1D_cState &U) {

  double rel_err, rel_err_last, precon, us(U[2]/U[1]);
  Levermore1D_cState U_temp(*this, us);
  Levermore1D_weights start(*this);
  ColumnVector A_step(Levermore1D_Vector::get_length());
  ColumnVector rhs(Levermore1D_Vector::get_length());
  ColumnVector Plambda(Levermore1D_Vector::get_length());
  DenseMatrix d2hda2(Levermore1D_Vector::get_length(),
			 Levermore1D_Vector::get_length());
  int count(0), count2(0);

  double damping(0.5); //0<=damping<1.0

  if(U_temp.relative_error(U) > 0.05) {
    cout << "!";
    MaxBoltz(U); //this is probably a better first guess if you're far away
    U_temp = Levermore1D_cState(*this,us);
  }

  for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
    rhs[i] = U[i+1]-U_temp[i+1];
  }
  rel_err = U_temp.relative_error(U);
  rel_err_last = rel_err;

  while( rel_err > 1e-10) { //fix tolerance later
    d2hda2 = U_temp.d2hda2(*this,us);

    if(count > 100) {
      cout << "============================================================"
	   << count << "   " << rel_err << endl 
	   << "U = " << U << endl
	   << "UA= " << U_temp << endl
	   << "A = " << (*this) << endl
	   << "Astart = " << start << endl
	   << d2hda2 << endl
	   << rhs << endl;
    }

    //precondition d2hda2 & rhs
    //for(int i=0;i<length;++i) {
    //  precon = d2hda2( ((i+1)/2)*2, 0 );
    //  precon = 1.0;
    //  for(int j=0;j<length;++j) {
    //	d2hda2(i,j) /= precon;
    //  }
    //  rhs(i) /= precon;
    //}
    Solve_LU_Decomposition(d2hda2,rhs,A_step);
    //A_step = rhs;
    //Solve_LS_Householder_F77(d2hda2,A_step,junk,length,length);
    //Solve_LAPACK_dgesv(d2hda2,A_step);
    if(count > 100) {
      cout << A_step << endl;
    }

    *this += A_step;
    U_temp.set_from_A(*this, us);
    rel_err = U_temp.relative_error(U);

    count2 = 0;
    while(rel_err >= rel_err_last) {
      cout << "@";
      A_step *= damping;
      *this -= A_step;
      U_temp.set_from_A(*this, us);
      rel_err = U_temp.relative_error(U);
      ++count2;
      if(count2 > 100) {
	cout << "Error, can't find step size to decrease rel_err" << endl;
	cout << "U = " << U << endl << "A = " << (*this) << endl 
	     << "UA= " << U_temp << endl << "Astart = " << start << endl
	     << "UAstart = " << Levermore1D_cState(start,us) << endl;
	MaxBoltz(U);
	U_temp = Levermore1D_cState(*this,us);
	rel_err = U_temp.relative_error(U);
	break;
      }
    }

    if(count == 50) {
      cout << "^";
      MaxBoltz(U); //try this?
      U_temp = Levermore1D_cState(*this,us);
      rel_err = U_temp.relative_error(U);
    }

    for(int i=0; i < Levermore1D_Vector::get_length(); ++i) {
      rhs[i] = U[i+1]-U_temp[i+1];
    }

    ++count;
    if(count > 200) return 1;
    rel_err_last = rel_err;
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
DenseMatrix Levermore1D_cState::d2hda2(const Levermore1D_weights &A, const double &us) const {

  double term;
  int j_min, j_max;
  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  for(int i = 0; i <= 2*(Levermore1D_Vector::get_length()-1); ++i) {
    term = moment(i,A,us);
    j_min = max(0,i-Levermore1D_Vector::get_length()+1);
    j_max = min(i,Levermore1D_Vector::get_length()-1);
    for(int j = j_min; j <= j_max; ++j)dm(i-j,j) = term;
  }

  return dm;
}

/********************************************************
 * Function: Levermore1D_cState::d2jda2                 *
 *                                                      *
 * Calculate Hessian of flux potential.                 *
 *                                                      *
 ********************************************************/
DenseMatrix Levermore1D_cState::d2jda2(const Levermore1D_weights &A, const double &us) const {

  double term;
  int j_min, j_max;
  DenseMatrix dm(Levermore1D_Vector::get_length(),
		 Levermore1D_Vector::get_length() );

  for(int i = 0; i <= 2*(Levermore1D_Vector::get_length()-1); ++i) {
    term = moment(i+1,A,us);
    j_min = max(0,i-Levermore1D_Vector::get_length()+1);
    j_max = min(i,Levermore1D_Vector::get_length()-1);
    for(int j = j_min; j <= j_max; ++j)dm(i-j,j) = term;
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
  double Wavel, Waver;

  Wavel = wavespeed_l*1.0;  //artificial damping
  Waver = wavespeed_r*1.0;

  if (Wavel >= ZERO) {
    Flux = Ul.F(Al);
  } else if (Waver <= ZERO) {
    Flux = Ur.F(Ar);
  } else {
    Flux =   ((Ul.F(Al)*Waver-Ur.F(Ar)*Wavel)
	      +(Ur-Ul)*(Wavel*Waver))/
             (Waver-Wavel);
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
			       const double &us_l,
			       const Levermore1D_weights &Ar,
			       const double &us_r) {
  Levermore1D_Vector Flux;

  for(int i=1; i<=Levermore1D_Vector::get_length(); ++i) {
    Flux[i] = Al.integrate_conserved_moment_pos(i,us_l);
    Flux[i] += Ar.integrate_conserved_moment_neg(i,us_r);
  }

  return Flux;
}
