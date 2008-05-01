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

  /* Constructors */
  BGK1D_Vector(void) : ColumnVector(BGK1D_Vector::get_length()) {}
  BGK1D_Vector(const BGK1D_Vector &B) : ColumnVector(B) {}
  BGK1D_Vector(const ColumnVector &C) : ColumnVector(C) {}

  /* Assignment operators */
  ColumnVector& operator=(const ColumnVector &V) {return ColumnVector::operator=(V);}
  ColumnVector& operator=(const BGK1D_Vector &V) {return ColumnVector::operator=(V);}

  /* Member Functions */
  void Maxwell_Boltzmann(double rho, double u, double v);
  ColumnVector velocities(void) const;
  double velocity(int n) const {return m_velocities[n];}
  double v_max(void) const;
  double v_min(void) const;

  /* Static Functions */
  static void setup(int l, double v_min, double v_max){
    double delta_v( (v_max-v_min)/(double)(l) );
    double v(v_min +(delta_v/2.0) );
    //ensure setup has not previously been done
    // and that v_max > v_min
    assert(m_set == 0 && v_max > v_min);
    m_set = 1;
    m_length = l;
    m_velocities = ColumnVector(l);
    for(int i=0; i<l; ++i) {
      m_velocities[i] = v;
      v += delta_v;
    }
  }

  static int get_length() {return m_length;}
  static int setup_done() {
    return m_set;
  }

 protected:

  static int m_length;
  static int m_set;
  static ColumnVector m_velocities;
};

#endif //_BGK1D_VECTOR_INCLUDED
