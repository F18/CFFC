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
  Levermore1D_Vector operator^(const Levermore1D_Vector &L) const;

  /* Functions */
  void copy_from(const Levermore1D_Vector &L);
  void zero() {set_all(0.0);}
  void one() {set_all(1.0);}
  void set_all(double in);
  void output(ostream &out) const;

  /* Static Functions */
  static void set_length(int l){
    assert(length_set == 0 || l == length);
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

inline Levermore1D_Vector Levermore1D_Vector::operator^(const Levermore1D_Vector &L) const{
  Levermore1D_Vector temp;
  for(int i=0;i<length;++i) {
    temp.m_values[i] = m_values[i] * L.m_values[i];
  }
  return temp;
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
  for(int i=0; i<length; ++i) {
    out << " " << m_values[i];
  }
}

/********************************************************
 *              External  Functions                     *
 ********************************************************/
inline ostream& operator<<(ostream &out, const Levermore1D_Vector &L) {
  L.output(out);
  return out;
}

#endif
