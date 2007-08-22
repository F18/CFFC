/* Third_order_tensor2D.h:  Header file defining third-order 2D Tensor class. */

#ifndef _THIRD_ORDER_TENSOR2D_INCLUDED
#define _THIRD_ORDER_TENSOR2D_INCLUDED

/* Include required C++ libraries. */

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cstdlib>

using namespace std;

/* Include math macro. */

#ifndef _MATH_MACROS_INCLUDED
#include "Math.h"
#endif // _MATH_MACROS_INCLUDED

/* Define the class. */

#define	NUM_COMP_TENSOR2D    3

/********************************************************
 * Class: Third_order_tensor2D (A symmetric tensor)     *
 *                                                      *
 * Member functions                                     *
 *      zero      -- Assign zero tensor.                *
 *      identity  -- Assign identity tensor.            *
 *      xxx       -- Return xxx-entry.                  *
 *      xxy       -- Return xxy-entry.                  *
 *      xyy       -- Return xyy-entry.                  *
 *      xxz       -- Return xxz-entry.                  *
 *      xzz       -- Return xzz-entry.                  *
 *      xyz       -- Return xyz-entry.                  *
 *      yyy       -- Return yyy-entry.                  *
 *      yyz       -- Return yyz-entry.                  *
 *      yzz       -- Return yzz-entry.                  *
 *      zzz       -- Return zzz-entry.                  *
 *                                                      *
 *  OPERATORS HAVE NOT YET BEEN IMPLEMENTED             *
 *  (safe to use default assignment operator)           *
 * Member operators                                     *
 *      T -- a third-order 2D tensor                    *
 *      a -- a scalar (double)                          *
 *                                                      *
 * T = T;                                               *
 * T = T + T;                                           *
 * T = T - T;                                           *
 * T = a * T;                                           *
 * T = T * a;                                           *
 * T = T / a;                                           *
 * T = +T;                                              *
 * T = -T;                                              *
 * T += T;                                              *
 * T -= T;                                              *
 * T += T;                                              *
 * T *= a;                                              *
 * T /= a;                                              *
 * T != T;                                              *
 * cout << T; (output function)                         *
 * cin  >> T; (input function)                          *
 *                                                      *
 ********************************************************/
class Third_order_tensor2D{
 private:
 public:
  double xxx,xxy,xyy,xzz,yyy,yzz;

    static double xxz,yyz,xyz,zzz; // Zero for two-D problems
    	                       // Made public so can access them.
			      
    /* Creation, copy, and assignment constructors. */
    Third_order_tensor2D(void) {
      zero();
    }

    Third_order_tensor2D(const Third_order_tensor2D &T) {
       xxx = T.xxx; xxy = T.xxy; xyy = T.xyy; xzz = T.xzz;
       yyy = T.yyy; yzz = T.yzz;
    }

    Third_order_tensor2D(const double &Txxx,
	     const double &Txxy,
	     const double &Txyy,
	     const double &Txzz,
	     const double &Tyyy,
	     const double &Tyzz) {
       xxx = Txxx; xxy = Txxy; xyy = Txyy; xzz = Txzz;
       yyy = Tyyy; yzz = Tyzz;
    }

    /* Destructor. */
    // ~Third_order_tensor2D(void);
    // Use automatically generated destructor.

    /* Zero (re-initialize) the tensor. */
    void zero(void);
    
    /* Input-output operators. */
    friend ostream &operator << (ostream &out_file, const Third_order_tensor2D &T);
    friend istream &operator >> (istream &in_file, Third_order_tensor2D &T);

};

/********************************************************
 * Third_order_tensor2D::zero -- Assign zero tensor.    *
 ********************************************************/
inline void Third_order_tensor2D::zero(void) {
    xxx = ZERO; xxy = ZERO; xyy = ZERO; xzz = ZERO;
    yyy = ZERO; yzz = ZERO;
}

/********************************************************
 * Third_order_tensor2D -- Input-output operators.      *
 ********************************************************/
inline ostream &operator << (ostream &out_file, const Third_order_tensor2D &T) {
  out_file.setf(ios::scientific);
  out_file << " " << T.xxx << " " << T.xxy << " " << T.xyy << " " << T.xzz
	   << " " << T.yyy << " " << T.yzz;
  out_file.unsetf(ios::scientific);
  return (out_file);
}

inline istream &operator >> (istream &in_file,Third_order_tensor2D &T) {
  in_file.setf(ios::skipws);
  in_file >> T.xxx >> T.xxy >> T.xyy >> T.xzz >> T.yyy >> T.yzz;
  in_file.unsetf(ios::skipws);
  return (in_file);
}

/********************************************************
 * Third_order_tensor2D -- Useful 2D tensor constants.  *
 ********************************************************/
const Third_order_tensor2D Third_order_tensor2D_ZERO(ZERO,ZERO,ZERO,ZERO,ZERO,ZERO);

#endif /* _TENSOR2D_INCLUDED  */
