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

/********************************************************
 * Class: Levermore1D_Vector                            *
 ********************************************************/
template<int N_moments>
class Levermore1D_Vector {
 protected:
  double m_values[N_moments];
 public:

  /* Default constructor */
  Levermore1D_Vector(void) {
    zero();
  }

  /* Copy constructor */
  Levermore1D_Vector(const Levermore1D_Vector &L) {
    copy_from(L);
  }

  /* Operators */

  double& operator[](int index);
  Levermore1D_Vector& operator=(const Levermore1D_Vector &L);
  Levermore1D_Vector operator+(const Levermore1D_Vector &L) const;
  Levermore1D_Vector& operator+=(const Levermore1D_Vector &L);
  Levermore1D_Vector operator-(const Levermore1D_Vector &L) const;
  Levermore1D_Vector& operator-=(const Levermore1D_Vector &L);
  double operator*(const Levermore1D_Vector &L) const;
  Levermore1D_Vector operator^(const Levermore1D_Vector &L) const;

  /* Functions */

  void copy_from(const Levermore1D_Vector &L);
  void zero();
  void one();
  void set_all(double in);
  void output(ostream &out) const;

};

/********************************************************
 *                  Operators                           *
 ********************************************************/
template<int N_moments>
inline double& Levermore1D_Vector<N_moments>::operator[](int index) {
  //I will start the index at (shudder) one to be consistent
  //with the rest of the code
  assert( index >= 1 && index <= N_moments );
  return m_values[index-1];
}

template<int N_moments>
inline Levermore1D_Vector<N_moments>& Levermore1D_Vector<N_moments>::operator=(const Levermore1D_Vector &L) {
  if(this == &L) {return (*this);}
  copy_from(L);
  return (*this);
}

template<int N_moments>
inline Levermore1D_Vector<N_moments> Levermore1D_Vector<N_moments>::operator+(const Levermore1D_Vector &L) const{
  Levermore1D_Vector<N_moments> temp;
  for(int i=0; i<N_moments; ++i) {
    temp.m_values[i] = m_values[i]+L.m_values[i];
  }
  return temp;
}

template<int N_moments>
inline Levermore1D_Vector<N_moments>& Levermore1D_Vector<N_moments>::operator+=(const Levermore1D_Vector &L) {
  for(int i=0; i<N_moments; ++i) {
    m_values[i] += L.m_values[i];
  }
  return (*this);
}

template<int N_moments>
inline Levermore1D_Vector<N_moments> Levermore1D_Vector<N_moments>::operator-(const Levermore1D_Vector &L) const{
  Levermore1D_Vector<N_moments> temp;
  for(int i=0; i<N_moments; ++i) {
    temp.m_values[i] = m_values[i]-L.m_values[i];
  }
  return temp;
}

template<int N_moments>
inline Levermore1D_Vector<N_moments>& Levermore1D_Vector<N_moments>::operator-=(const Levermore1D_Vector &L) {
  for(int i=0; i<N_moments; ++i) {
    m_values[i] -= L.m_values[i];
  }
  return (*this);
}

template<int N_moments>
inline double Levermore1D_Vector<N_moments>::operator*(const Levermore1D_Vector &L) const{
  double temp;
  for(int i=0; i<N_moments; ++i) {
    temp += m_values[i]*L.m_values[i];
  }
  return temp;
}

template<int N_moments>
inline Levermore1D_Vector<N_moments> Levermore1D_Vector<N_moments>::operator^(const Levermore1D_Vector &L) const{
  Levermore1D_Vector<N_moments> temp;
  for(int i=0; i<N_moments; ++i) {
    temp.m_values[i] = m_values[i]*L.m_values[i];
  }
  return temp;
}

/********************************************************
 *                  Functions                           *
 ********************************************************/
template<int N_moments>
inline void Levermore1D_Vector<N_moments>::copy_from(const Levermore1D_Vector &L) {
  for(int i=0; i<N_moments; ++i) {
    m_values[i] = L.m_values[i];
  }
}

template<int N_moments>
inline void Levermore1D_Vector<N_moments>::zero() {
  set_all(0.0);
}

template<int N_moments>
inline void Levermore1D_Vector<N_moments>::one() {
  set_all(1.0);
}

template<int N_moments>
inline void Levermore1D_Vector<N_moments>::set_all(double in) {
  for(int i=0; i<N_moments; ++i) {
    m_values[i] = in;
  }
}

template<int N_moments>
inline void Levermore1D_Vector<N_moments>::output(ostream &out) const{
  for(int i=0; i<N_moments; ++i) {
    out << " " << m_values[i];
  }
}

/********************************************************
 *              External  Functions                     *
 ********************************************************/
template<int N_moments>
inline ostream& operator<<(ostream &out, const Levermore1D_Vector<N_moments> &L) {
  L.output(out);
  return out;
}

#endif
