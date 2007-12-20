#ifndef _ARRAY1D_MATH_INCLUDED
#define _ARRAY1D_MATH_INCLUDED


/*==============================================================================
 * Included Files
 *============================================================================*/

//----- Standard Library -----//

#include<cstddef>


namespace Math
{

namespace Array1D
{


/*******************************************************************************
 *
 * struct Vector_Subscript
 *
 * Purpose
 * =======
 *
 *   Implments a structure for indexing a dimension of an array.  This contains
 *   a starting index, an ending index, and a stride.
 *
 ******************************************************************************/

struct Vector_Subscript

{


/*==============================================================================
 * Data members
 *============================================================================*/

  public:

   int                      min;  // First subscript
   int                      max;  // Last subscript
   int                   stride;  // Stride


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

   public:

//--Default constructor

   Vector_Subscript() { }

//--Constructor

   Vector_Subscript(const int imin, const int imax, const int istride):
      min(imin), max(imax), stride(istride) { }

//--Copy constructor

   Vector_Subscript(const Vector_Subscript& vs):
      min(vs.min), max(vs.max), stride(vs.stride) { }

//--Assignment operator

   Vector_Subscript& operator=(const Vector_Subscript& vs) {
      if (&vs != this) {
         min = vs.min;
         max = vs.max;
         stride = vs.stride;
      }
      return *this;
   }

//--Use automatic destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return number of elements in the subscript

   int nelem() const {
      return (max - min + 1)/stride;
   }

//--Return exit mark for loop termination (e.g. i != vs.exit())

   int exit() const {
      return max + stride;
   }

};


/*******************************************************************************
 *
 * Functions on 1D arrays.
 *
 ******************************************************************************/


/*==============================================================================
 *
 * function maxval
 *
 * Purpose
 * =======
 *
 *   Finds the maximum value in an array of size n (first if more than 1)
 *
 * Notes
 * =====
 *
 *   > must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline const T& maxval(const T* p, int n)

{

   const T* q = p;
   while ( --n ) if ( *(++p) > *q ) q = p;
   return *q;

}


/*==============================================================================
 *
 * function maxval
 *
 * Purpose
 * =======
 *
 *   Finds the maximum value in an array of size n (first if more than 1).  T
 *   is usually a user-defined type and Cmp is a functor that returns a special
 *   less than comparison bewteen two T's x and y when operator(x, y) is called.
 *
 * Notes
 * =====
 *
 *   n > 0
 *
 *============================================================================*/

template <typename T, typename Cmp>
inline const T& maxval(const T* p, int n, const Cmp& cmp)

{

   const T* q = p;
   while ( --n ) if ( cmp(*q, *(++p)) ) q = p;
   return *q;

}


/*==============================================================================
 *
 * function minval
 *
 * Purpose
 * =======
 *
 *   Finds the minimum value in an array of size n (first if more than 1)
 *
 * Notes
 * =====
 *
 *   < must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline const T& minval(const T* p, size_t n)

{

    const T* q = p;
    while ( --n ) if ( *(++p) < *q ) q = p;
    return *q;

}


/*==============================================================================
 *
 * function minval
 *
 * Purpose
 * =======
 *
 *   Finds the minimum value in an array of size n (first if more than 1).  T
 *   is usually a user-defined type and Cmp is a functor that returns a special
 *   less than comparison bewteen two T's x and y when operator(x, y) is called.
 *
 * Notes
 * =====
 *
 *   n > 0
 *
 *============================================================================*/

template <typename T, typename Cmp>
inline const T& minval(const T* p, size_t n, const Cmp& cmp)

{

   const T* q = p;
   while ( --n ) if ( cmp(*(++p), *q) ) q = p;
   return *q;

}


/*==============================================================================
 *
 * function maxloc
 *
 * Purpose
 * =======
 *
 *   Finds the location of the maximum value in an array of size n (first if
 *   more than 1)
 *
 * Notes
 * =====
 *
 *   < must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline size_t maxloc(const T* p, const size_t N)

{

   size_t i = N;
   const T* q = p;
   for ( size_t n = N ; --n ; ) {
      if ( *q < *(++p) ) {
         q = p;
         i = n;
      }
   }
   return N - i;

}


/*==============================================================================
 *
 * function maxloc
 *
 * Purpose
 * =======
 *
 *   Finds the location of the maximum value in an array of size n (first if
 *   more than 1).  T is usually a user-defined type and Cmp is a functor that
 *   returns a special less than comparison bewteen two T's x and y when
 *   operator(x, y) is called.
 *
 * Notes
 * =====
 *
 *   n > 0
 *
 *============================================================================*/

template <typename T, typename Cmp>
inline size_t maxloc(const T* p, const size_t N, const Cmp& cmp)

{

   size_t i = N;
   const T* q = p;
   for ( size_t n = N ; --n ; ) {
      if ( cmp(*q, *(++p)) ) {
         q = p;
         i = n;
      }
   }
   return N - i;

}


/*==============================================================================
 *
 * function minloc
 *
 * Purpose
 * =======
 *
 *   Finds the location of the minimum value in an array of size n (first if
 *   more than 1)
 *
 * Notes
 * =====
 *
 *   < must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline size_t minloc(const T* p, const size_t N)

{

   size_t i = N;
   const T* q = p;
   for ( size_t n = N ; --n ; ) {
      if ( *(++p) < *q ) {
         q = p;
         i = n;
      }
   }
   return N - i;

}


/*==============================================================================
 *
 * function minloc
 *
 * Purpose
 * =======
 *
 *   Finds the location of the minimum value in an array of size n (first if
 *   more than 1).  T is usually a user-defined type and Cmp is a functor that
 *   returns a special less than comparison bewteen two T's x and y when
 *   operator(x, y) is called.
 *
 * Notes
 * =====
 *
 *   n > 0
 *
 *============================================================================*/

template <typename T, typename Cmp>
inline size_t minloc(const T* p, const size_t N, const Cmp& cmp)

{

   size_t i = N;
   const T* q = p;
   for ( size_t n = N ; --n ; ) {
      if ( cmp(*(++p), *q) ) {
         q = p;
         i = n;
      }
   }
   return N - i;

}


/*==============================================================================
 *
 * function sum
 *
 * Purpose
 * =======
 *
 *   Sums n elements in the array
 *
 * Notes
 * =====
 *
 *   += must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline T sum(const T* p, int n)

{

    T total = *p;
    while ( --n ) total += *(++p);
    return total;

}


/*==============================================================================
 *
 * function product
 *
 * Purpose
 * =======
 *
 *   Mutliplies n elements in the array
 *
 * Notes
 * =====
 *
 *   *= must be defined for T
 *   n > 0
 *
 *============================================================================*/

template <typename T>
inline T product(const T* p, int n)

{

    T total = *p;
    while ( --n ) total *= *(++p);
    return total;

}


/*==============================================================================
 *
 * function accumulate
 *
 * Purpose
 * =======
 *
 *   Accumulates the objects in an array using functor Accum.  Accum must
 *   return a result of type T when operator (x, y) is invoked on two Ts.  Most
 *   commonly, Accum will sum or multiply the objects of the array
 *
 * Notes
 * =====
 *
 *   n > 0
 *
 *============================================================================*/

template <typename T, typename Accum>
inline T accumulate(const T* p, int n, const Accum& accum)

{

    T total = *p;
    while ( --n ) total = accum(total, *(++p));
    return total;

}


}  // Close namespace Array1D

}  // Close namespace Math

#endif // _ARRAY1D_MATH_INCLUDED
