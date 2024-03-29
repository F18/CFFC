/* Levermore1DVector.h:  Header file defining 1D Levermore Vector Class. */

#ifndef _LEVERMORE1D_VECTOR_INCLUDED
#define _LEVERMORE1D_VECTOR_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>

using namespace std;

/* Include Headers */
#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif //_MATRIX_INCLUDED

/********************************************************
 * Class: Levermore1D_Vector                            *
 ********************************************************/
class Levermore1D_Vector {
 protected:
  double *m_values;
  static int length;
  static int length_set;
 public:

  /* Default constructor */
  Levermore1D_Vector(void) {
    m_values = new double[length];
  }

  /* Copy constructor */
  Levermore1D_Vector(const Levermore1D_Vector &L) {
    m_values = new double[length];
    copy_from(L);
  }

  /* Constructor from column vector */
  Levermore1D_Vector(const ColumnVector &CV) {
    assert(length == CV.size());
    m_values = new double[length];
    for(int i=0;i<length;++i) {
      m_values[i] = CV[i];
    }
  }

  /* Destructor */
  ~Levermore1D_Vector() {delete [] m_values;}

  /* Operators */
  double& operator[](int index);
  const double& operator[](int index) const;
  Levermore1D_Vector& operator=(const Levermore1D_Vector &L);
  Levermore1D_Vector operator+(const Levermore1D_Vector &L) const;
  Levermore1D_Vector& operator+=(const Levermore1D_Vector &L);
  Levermore1D_Vector operator-(const Levermore1D_Vector &L) const;
  Levermore1D_Vector& operator-=(const Levermore1D_Vector &L);
  double operator*(const Levermore1D_Vector &L) const;
  Levermore1D_Vector operator*(const double &d) const;
  Levermore1D_Vector operator/(const double &d) const;
  Levermore1D_Vector operator^(const Levermore1D_Vector &L) const;

  /* Functions */
  void copy_from(const Levermore1D_Vector &L);
  void zero() {set_all(0.0);}
  void one() {set_all(1.0);}
  void set_all(double in);
  void output(ostream &out) const;
  void input(istream &in) const;

  /* Static Functions */
  static void set_length(int l){
    //ensure length has not previously been set
    assert(length_set == 0 || l == length);
    //ensure this is a valid length
    assert(l > 2 && l%2 ==1);

    length = l;
    length_set = 1;
  }

  static int get_length() {return length;}

  static int length_is_set() {
    return length_set;
  }

};

/********************************************************
 *                  Operators                           *
 ********************************************************/
inline double& Levermore1D_Vector::operator[](int index) {
  //I will start the index at (shudder) one to be consistent
  //with the rest of the code
  assert( index >= 1 && index <= length );
  return m_values[index-1];
}

inline const double& Levermore1D_Vector::operator[](int index) const {
  //I will start the index at (shudder) one to be consistent
  //with the rest of the code
  assert( index >= 1 && index <= length );
  return m_values[index-1];
}

inline Levermore1D_Vector& Levermore1D_Vector::operator=(const Levermore1D_Vector &L) {
  //if(this == &L) {return (*this);}
  copy_from(L);
  return (*this);
}

inline Levermore1D_Vector Levermore1D_Vector::operator+(const Levermore1D_Vector &L) const {
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] + L.m_values[i];
  }
  return temp;
}

inline Levermore1D_Vector& Levermore1D_Vector::operator+=(const Levermore1D_Vector &L) {
  for(int i=0;i<length;++i) {
    m_values[i] += L.m_values[i];
  }
  return (*this);
}

inline Levermore1D_Vector Levermore1D_Vector::operator-(const Levermore1D_Vector &L) const{
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] - L.m_values[i];
  }
  return temp;
}

inline Levermore1D_Vector& Levermore1D_Vector::operator-=(const Levermore1D_Vector &L) {
  for(int i=0;i<length;++i) {
    m_values[i] -= L.m_values[i];
  }
  return (*this);
}

inline double Levermore1D_Vector::operator*(const Levermore1D_Vector &L) const{
  double temp(0.0);
  for(int i=0;i<length;++i) {
    temp += m_values[i] * L.m_values[i];
  }
  return temp;
}

inline Levermore1D_Vector Levermore1D_Vector::operator*(const double &d) const{
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] * d;
  }
  return temp;
}

inline Levermore1D_Vector Levermore1D_Vector::operator/(const double &d) const{
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] / d;
  }
  return temp;
}

inline Levermore1D_Vector Levermore1D_Vector::operator^(const Levermore1D_Vector &L) const{
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] * L.m_values[i];
  }
  return temp;
}

/********************************************************
 *                 External Operators                   *
 ********************************************************/
inline Levermore1D_Vector operator *(const DenseMatrix &M, const Levermore1D_Vector &LVec) {
  assert(M.size(0) ==  LVec.get_length() &&
	 M.size(0) == M.size(1));
  int i, j; Levermore1D_Vector lv;
  lv.zero();
  for ( i = 0; i <= M.size(0)-1; ++i ) {
    for ( j = 0 ; j <= M.size(1)-1; ++j ) lv[i+1]+=M(i,j)*LVec[j+1];
  } /* endfor */
  return (lv);
}


/********************************************************
 *                  Functions                           *
 ********************************************************/
inline void Levermore1D_Vector::copy_from(const Levermore1D_Vector &L) {
  for(int i=0;i<length;++i) {
    m_values[i] = L.m_values[i];
  }
}

inline void Levermore1D_Vector::set_all(double in) {
  for(int i=0;i<length;++i) {
    m_values[i] = in;
  }
}

inline void Levermore1D_Vector::output(ostream &out) const{
  out.precision(16);
  for(int i=0; i<length; ++i) {
    out << " " << m_values[i];
  }
  out.precision(6);
}

inline void Levermore1D_Vector::input(istream &in) const{
  for(int i=0; i<length; ++i) {
    in >> m_values[i];
  }
}

/********************************************************
 *              External  Functions                     *
 ********************************************************/
inline ostream& operator<<(ostream &out, const Levermore1D_Vector &L) {
  L.output(out);
  return out;
}

inline istream& operator>>(istream &in, const Levermore1D_Vector &L) {
  L.input(in);
  return in;
}

#endif
