/* BGK1D_Vector.h:  Header file for vector defining discretization of 
                    distribution function for BGK1D solver. */

#ifndef _BGK1D_VECTOR_INCLUDED
#define _BGK1D_VECTOR_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <cmath>

using namespace std;

#ifndef _MATRIX_INCLUDED
#include "../Math/Matrix.h"
#endif //_MATRIX_INCLUDED

#ifndef _LINEARSYSTEMS_INCLUDED
#include "../Math/LinearSystems.h"
#endif //_LINEARSYSTEMS_INCLUDED

#ifndef _GAS_CONSTANTS_INCLUDED
#include "../Physics/GasConstants.h"
#endif // _GAS_CONSTANTS_INCLUDED

class BGK1D_Vector : public ColumnVector {

 public:

  BGK1D_Vector(void) : ColumnVector(BGK1D_Vector::get_length()) {}
  BGK1D_Vector(const BGK1D_Vector &B) : ColumnVector(B) {}
  BGK1D_Vector(const ColumnVector &C) : ColumnVector(C) {}

  /* Assignment operators */
  ColumnVector& operator=(const ColumnVector &V) {return ColumnVector::operator=(V);}
  ColumnVector& operator=(const BGK1D_Vector &V) {return ColumnVector::operator=(V);}

  /* Static Functions */
  static void set_length(int l){
    //ensure length has not previously been set
    assert(m_length_set == 0 || l == m_length);
    m_length = l;
    m_length_set = 1;
  }

  static int get_length() {return m_length;}

  static int length_is_set() {
    return m_length_set;
  }


 protected:

  static int m_length;
  static int m_length_set;
};

#endif //_BGK1D_VECTOR_INCLUDED
