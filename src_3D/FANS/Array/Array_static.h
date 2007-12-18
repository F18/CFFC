#ifndef _STATIC_ARRAY_INCLUDED
#define _STATIC_ARRAY_INCLUDED


/*==============================================================================
 * Included Files
 *============================================================================*/

//----- Standard Library -----//

#include<iostream>
#include<stdexcept>
#include<string>

//===== End of Includes =====//

/*==============================================================================
 * Forward declarations of math matrix types.  Required for declarations of
 * friends (see Static with rank 2).
 *============================================================================*/

namespace Math 
{
 
   template<typename Vec, typename Vec2, size_t Dim, size_t N,
      Array::ArraySizeTr_t ArraySizeTr>
      struct AssignRepPolicy;
   // Binary expressions
   template <typename T1, typename T2, typename OP1, typename OP2>
      class Array_Add;
   template <typename T1, typename T2, typename OP1, typename OP2>
      class Array_Subtract;
   template <typename T1, typename T2, typename OP1, typename OP2>
      class Array_Mult;
   template <typename T1, typename T2, typename OP1, typename OP2>
      class Array_Div;
   // Unary expressions
   template <typename MT, typename OP1>
      class Array_Neg;
   template <typename MT, typename OP1>
      class Array_Abs;

}


namespace Array
{


/*==============================================================================
 * Forward declarations of sister array types
 *============================================================================*/

//--Dynamic

template <typename T, int Rank>
class Dynamic;
template <typename T>
class Dynamic_Base;

//--Reference

template <typename T, int Rank>
class Reference;
template <typename T>
class Reference_Base;


/*******************************************************************************
 *
 * class Static
 *
 * Purpose
 * =======
 *
 *   Implements a static array (an array for which the dimensions are known
 *   at compiler time) with up to 7 dimensions by indexing a 1 dimensional
 *   array.
 *
 * Constructors
 * ============
 *
 *   The class takes from 2 to 8 template arguments.  The first is the type of
 *   data in the array and the remaining are the dimensions of the array.
 *
 *   Static()          -- default constructor
 *
 *   Static(const Static& A)
 *                     -- copy constructor for argument of same type, rank, and
 *                        dimensions
 *     A                  (I) array to copy
 *
 *   Static& operator=(const Static& A)
 *                     -- assignment operator for argument of same type, rank,
 *                        and dimensions
 *     A                  (I) array to copy
 *     return             (O) updated reference
 *
 *   Static& operator=(const Static<T2>& A)
 *                     -- assignment operator for argument of same rank and
 *                        dimensions but different type. 
 *     A                  (I) array to copy
 *     return             (O) updated reference
 *
 * Member Functions
 * ================
 *
 *   const T& max()    -- max value in the array ( < must be defined for T)
 *     return             (O) Reference to maximum value in array
 *
 *   const T& max(cmp) -- max value using predicate cmp
 *     cmp                (I) Predicate defining < for two elements
 *     return             (O) Reference to maximum value in array
 *
 *   const T& min()    -- min value in the array ( < must be defined for T)
 *     return             (O) Reference to minimum value in the array
 *
 *   const T& min(cmp) -- min value using predicate cmp
 *     cmp                (I) Predicate defining < for two elements
 *     return             (O) Reference to minimum value in array
 *
 *   T sum()           -- sum of all values in array ( += must be defined for T)
 *     return             (O) Sum
 *
 *   T product()       -- product of all values in array ( *= must be defined
 *                        for T)
 *     return             (O) Product
 *
 *   T accumulate(accum)
 *                     -- Accumulate using functor accum
 *     accum              (I) Functor which accumulates two elements and returns
 *                            result
 *     return             (O) Total accumulation
 *
 *   size_t dim_size(const size_t dim = 0)
 *                     -- returns length of requested dimension or total size
 *                        if dim = 0.  If that dimension does not exists, 0 is
 *                        returned
 *     dim                (I) dimension to return size of.  A value of 0 will
 *                            total size of the array
 *     return             (O) size of the dimension
 *
 * Member Operators
 * ================
 *
 *   T& operator()(const size_t i0 ... i6)
 *                     -- indexes the array (number if indexes = rank of array)
 *     i0 ... i6          (I) indices to the array (number of indices = rank)
 *     return             (O) reference to array element
 *
 *   operator=(const T2& val)
 *                     -- assigns "val" to each index in the array (note that
 *                        "val" may be of a different type)
 *     val                (I) value to assign
 *
 *   ostream& operator<<(ostream& os, const Static<T, ...>& A)
 *                     -- writes out the array elements in storage order with
 *                        single space between each element (assumes << defined
 *                        for T)
 *     os                 (I) output stream
 *     A                  (I) static array
 *     return             (O) output stream
 *
 * Data Members
 * ============
 *
 *   size              -- Overall size of the array
 *   data              -- Pointer to the data
 *
 * Notes
 * =====
 *
 *   - The dimensions of the array must be known at compiler time
 *   - Uses C++ storage order
 *   - Define _STATIC_ARRAY_CHECK_BOUNDS to enable array bounds checking.  If
 *     a bound is exceeded, the value of the incorrect index is printed and the
 *     program aborts.
 *   - The rank 7 class is the general class and all other ranks are
 *     specializations.
 *   - The class uses template metaprogramming to unroll loops for arrays that
 *     have a size less than or equal to 60.  This behaviour can be modified in
 *     file "Array_traits.h"
 *
 * Exceptions
 * ==========
 *
 *   std::out_of_range
 *
 * Programmer
 * ==========
 *
 *   S. Guzik         13 October  2004    - Original code
 *   S. Guzik         09 November 2005    - Modified to use template
 *                                        metaprogramming
 *                                        - Restructured with base class
 *
 * Example
 * =======
 *
 *------------------------------------------------------------------------------

#include<iostream>
#include "static_array.h"

int main() {

   Array::Static<double, 3, 4> a;
   a = 2.;
   a(2, 2) = 4.;
   std::cout << "sum: " << a.sum() << std::endl;
   std::cout << "max: " << a.max() << std::endl;
   std::cout << "min: " << a.min() << std::endl;
   Array::Static<double, 3, 4> b = a;
   std::cout << "b[2][1]: " << b(2, 1) << " b[2][2]: " << b(2, 2) << std::endl;
   std::cout << "All of b: " << b << std::endl;

}

 *------------------------------------------------------------------------------
 *
 ******************************************************************************/


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static_Base
 *
 * Purpose
 * =======
 *
 *   Provides some base data and functions used by all of the Dynamic
 *   classes.
 *
 ******************************************************************************/

template<class T, size_t Size>
class Static_Base
   
{


/*==============================================================================
 * Public types
 *============================================================================*/

  public:

   typedef T value_type;


/*==============================================================================
 * Public member functions
 *============================================================================*/

  public:

//--Max

   const T& max() const {
      return MaxValPolicy<T, Size, Traits<Size>::size, std::less<T> >::eval(
         _data, std::less<T>());
   }

   // With user defined less than predicate
   template<typename Cmp>
   const T& max(const Cmp& cmp) const {
      return MaxValPolicy<T, Size, Traits<Size>::size, Cmp>::eval(_data, cmp);
   }

//--Min

   const T& min() const {
      return MinValPolicy<T, Size, Traits<Size>::size, std::less<T> >::eval(
         _data, std::less<T>());
   }

   // With user defined less than predicate
   template<typename Cmp>
   const T& min(const Cmp& cmp) const {
      return MinValPolicy<T, Size, Traits<Size>::size, Cmp>::eval(_data, cmp);
   }

//--Sum

   T sum() const {
      return AccumulatePolicy<T, Size, Traits<Size>::size, std::plus<T> >::eval(
         _data, std::plus<T>());
   };

//--Product

   T product() const {
      return AccumulatePolicy<T, Size, Traits<Size>::size,
         std::multiplies<T> >::eval(_data, std::multiplies<T>());
   };

//--Accumulate using user defined accumulation functor

   template<typename Accum>
   T accumulate(const Accum& accum) const {
      return AccumulatePolicy<T, Size, Traits<Size>::size, Accum>::eval(
         _data, accum);
   }

//--Return pointer to beginning of array

   T *begin() { return _data; }
   const T *begin() const { return _data; }

//--Return pointer to end of array

   T *end() { return _data + Size; }
   const T *end() const { return _data + Size; }

//--Return size of the array

   size_t size() const { return Size; }


/*==============================================================================
 * Protected constructors and destructors
 *============================================================================*/

  protected:

//--Default constructor

   Static_Base() { }

//--Copy constructor

   Static_Base(const Static_Base& A) { assign_base(A); }

//--Destructor

   ~Static_Base() { }


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assignment of base
 *--------------------------------------------------------------------*/

   // Same type
   void assign_base(const Static_Base &A)
   {
      AssignPolicy<T, T, Size, Traits<Size>::size>::eval(_data, A._data);
   }

   // Different types and/or size
   template<typename T2, size_t Size2>
   void assign_base(const Static_Base<T2, Size2> &A)
   {
      enum { MinSize = (Size < Size2) ? Size : Size2 };
      AssignPolicy<T, T2, MinSize, Traits<MinSize>::size>::eval(_data, A._data);
   }

   // From a dynamic array
   template <typename T2>
   void assign_base(const Dynamic_Base<T2> &A)
   {
      T *p = _data;
      const T2 *a = A._data;
      for ( size_t n = std::min(Size, A.size) ; n-- ; ) *p++ = *a++;
   }

   // From a reference array
   template <typename T2>
   void assign_base(const Reference_Base<T2> &A)
   {
      T *p = _data;
      const T2 *a = A._data;
      for ( size_t n = std::min(Size, A.size) ; n-- ; ) *p++ = *a++;
   }

/*--------------------------------------------------------------------*
 * Assignment of constant to base
 *--------------------------------------------------------------------*/

   template<typename T2>
   void assign_constant(const T2 &val)
   {
      AssignConstantPolicy<T, T2, Size, Traits<Size>::size>::eval(_data, val);
   }


/*==============================================================================
 * Private constructors
 *============================================================================*/

  private:

   Static_Base& operator=(const Static_Base&);
                                        // Assignment not allowed


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:

   T _data[Size];                       // Data for the array


/*==============================================================================
 * Tools for debugging
 *============================================================================*/

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
   void check_bounds(const int dim, const size_t index, const size_t iMax) const
   {
      if ( index < 0 || index >= iMax ) {
         char oText[128];
         std::sprintf(oText, "Index %d is outside array bounds:\n"
                      "Index value: %d\n"
                      "Valid range: 0 : %d\n"
                      "Error in class Array::Static", dim, index, iMax-1);
         throw std::out_of_range(oText);
      }
   }
#endif


/*==============================================================================
 * Support for MPI
 *============================================================================*/

#ifdef USE_MPI
  public:

   // Type T is assumed to store without gaps in an array.  If there are gaps,
   // the gaps must be represented in type MPIuT.
   MPI::Datatype mpi_struct(const MPI::Datatype &MPIuT) const
   {
      return MPIuT.Create_contiguous(Size);
   }

   void *mpi_ref()
   {
      return _data;
   }

   const void *mpi_ref() const
   {
      return _data;
   }
#endif

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 7
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 7.  This is the general template class
 *   and specializations are used for arrays of rank 6 - 1.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1 = 0, size_t dim2 = 0,
         size_t dim3 = 0, size_t dim4 = 0, size_t dim5 = 0, size_t dim6 = 0>
class Static: public Static_Base<T, dim0*dim1*dim2*dim3*dim4*dim5*dim6>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() :
      Static_Base<T, _size>(),
      stride5(dim6),
      stride4(stride5*dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Copy constructor

   Static(const Static& A) :
      Static_Base<T, _size>(A),
      stride5(dim6),
      stride4(stride5*dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      case 3:
         return dim2;
         break;
      case 4:
         return dim3;
         break;
      case 5:
         return dim4;
         break;
      case 6:
         return dim5;
         break;
      case 7:
         return dim6;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4, const size_t i5,
                 const size_t i6)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
      this->check_bounds(5, i5, dim5);
      this->check_bounds(6, i6, dim6);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5*stride5 + i6];
   }

   const T& operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4, const size_t i5,
                       const size_t i6) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
      this->check_bounds(5, i5, dim5);
      this->check_bounds(6, i6, dim6);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5*stride5 + i6];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride5;                // Stride for index 5
   const size_t stride4;                // Stride for index 4
   const size_t stride3;                // Stride for index 3
   const size_t stride2;                // Stride for index 2
   const size_t stride1;                // Stride for index 1
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1*dim2*dim3*dim4*dim5*dim6;
                                        // Total size of the array


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 7>;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 6
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 6.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1, size_t dim2 , size_t dim3,
         size_t dim4, size_t dim5>
class Static<T, dim0, dim1, dim2, dim3, dim4, dim5, 0>:
   public Static_Base<T, dim0*dim1*dim2*dim3*dim4*dim5>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() :
      Static_Base<T, _size>(),
      stride4(dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Copy constructor

   Static(const Static& A) :
      Static_Base<T, _size>(A),
      stride4(dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      case 3:
         return dim2;
         break;
      case 4:
         return dim3;
         break;
      case 5:
         return dim4;
         break;
      case 6:
         return dim5;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4, const size_t i5)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
      this->check_bounds(5, i5, dim5);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5];
   }

   const T& operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4, const size_t i5) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
      this->check_bounds(5, i5, dim5);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride4;                // Stride for index 4
   const size_t stride3;                // Stride for index 3
   const size_t stride2;                // Stride for index 2
   const size_t stride1;                // Stride for index 1
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1*dim2*dim3*dim4*dim5;
                                        // Total size of the array


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 6>;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 5
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 5.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1, size_t dim2 , size_t dim3,
         size_t dim4>
class Static<T, dim0, dim1, dim2, dim3, dim4, 0, 0>:
   public Static_Base<T, dim0*dim1*dim2*dim3*dim4>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() :
      Static_Base<T, _size>(),
      stride3(dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Copy constructor

   Static(const Static& A) :
      Static_Base<T, _size>(A),
      stride3(dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      case 3:
         return dim2;
         break;
      case 4:
         return dim3;
         break;
      case 5:
         return dim4;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4];
   }

   const T& operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
      this->check_bounds(4, i4, dim4);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride3;                // Stride for index 3
   const size_t stride2;                // Stride for index 2
   const size_t stride1;                // Stride for index 1
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1*dim2*dim3*dim4;
                                        // Total size of the array


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 5>;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 4
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 4.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1, size_t dim2 , size_t dim3>
class Static<T, dim0, dim1, dim2, dim3, 0, 0, 0>:
   public Static_Base<T, dim0*dim1*dim2*dim3>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() :
      Static_Base<T, _size>(),
      stride2(dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Copy constructor

   Static(const Static& A) :
      Static_Base<T, _size>(A),
      stride2(dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1) {
   }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      case 3:
         return dim2;
         break;
      case 4:
         return dim3;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3];
   }

   const T& operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
      this->check_bounds(3, i3, dim3);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride2;                // Stride for index 2
   const size_t stride1;                // Stride for index 1
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1*dim2*dim3;
                                        // Total size of the array


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 4>;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 3
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 3.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1, size_t dim2>
class Static<T, dim0, dim1, dim2, 0, 0, 0, 0>:
   public Static_Base<T, dim0*dim1*dim2>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() :
      Static_Base<T, _size>(),
      stride1(dim2),
      stride0(stride1*dim1) {
   }

//--Copy constructor

   Static(const Static& A) :
      Static_Base<T, _size>(A),
      stride1(dim2),
      stride0(stride1*dim1) {
   }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      case 3:
         return dim2;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1, const size_t i2)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2];
   }

   const T& operator()(const size_t i0, const size_t i1, const size_t i2) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
      this->check_bounds(2, i2, dim2);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride1;                // Stride for index 1
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1*dim2;
                                        // Total size of the array


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 3>;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 2
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 2.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0, size_t dim1>
class Static<T, dim0, dim1, 0, 0, 0, 0, 0>: public Static_Base<T, dim0*dim1>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() : Static_Base<T, _size>(), stride0(dim1) { }

//--Copy constructor

   Static(const Static& A) : Static_Base<T, _size>(A), stride0(dim1) { }

//--Assignment operator

   // Same type and size
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use sythesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
         return _size;
         break;
      case 1:
         return dim0;
         break;
      case 2:
         return dim1;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0, const size_t i1)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
#endif

      return this->_data[i0*stride0 + i1];
   }

   const T& operator()(const size_t i0, const size_t i1) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
      this->check_bounds(1, i1, dim1);
#endif

      return this->_data[i0*stride0 + i1];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   const size_t stride0;                // Stride for index 0
   static const size_t _size = dim0*dim1;
                                        // Total size of the array


/*==============================================================================
 * Private member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * This single indexing operator is for the expression templates
 * defined for matrices in the Math classes.  This should be kept
 * private because it can lead to very difficult-to-find bugs in user
 * code otherwise.  Note, however, that an extensive list of friends
 * is required, basically all the math matrix expressions.
 *--------------------------------------------------------------------*/

  private:
   T &operator()(const size_t i) { return this->_data[i]; }
   const T &operator()(const size_t i) const { return this->_data[i]; }


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 2>;

/*--------------------------------------------------------------------*
 * Friends for the math matrix assignment policy and for all math
 * matrix expressions.
 *--------------------------------------------------------------------*/

   template<typename Vec, typename Vec2, size_t Dim, size_t N,
      ArraySizeTr_t ArraySizeTr>
      friend struct Math::AssignRepPolicy;
   // Binary expressions
   template <typename T1, typename T2, typename OP1, typename OP2>
      friend class Math::Array_Add;
   template <typename T1, typename T2, typename OP1, typename OP2>
      friend class Math::Array_Subtract;
   template <typename T1, typename T2, typename OP1, typename OP2>
      friend class Math::Array_Mult;
   template <typename T1, typename T2, typename OP1, typename OP2>
      friend class Math::Array_Div;
   // Unary expressions
   template <typename MT, typename OP1>
      friend class Math::Array_Neg;
   template <typename MT, typename OP1>
      friend class Math::Array_Abs;

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Static - rank 1
 *
 * Purpose
 * =======
 *
 *   Implements a static array of rank 1.  This is a specialized template class
 *   of the general rank 7 template class.
 *
 ******************************************************************************/

template<typename T, size_t dim0>
class Static<T, dim0, 0, 0, 0, 0, 0, 0>: public Static_Base<T, dim0>

{

  public:


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

//--Default constructor

   Static() : Static_Base<T, _size>() { }

//--Copy constructor

   Static(const Static& A) : Static_Base<T, _size>(A) { }

//--Assignment operator

   // Same type
   Static& operator=(const Static& A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or size
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Static& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                     SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a dynamic array
   template<typename T2, int Rank>
   Static& operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Static& operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

//--Use synthesized destructor


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Return size of dimensions

   size_t dim_size(const size_t dim = 0) {
      switch ( dim ) {
      case 0:
      case 1:
         return dim0;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Public operators
 *============================================================================*/

//--Index the array

   T& operator()(const size_t i0)
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
#endif

      return this->_data[i0];
   }

   const T& operator()(const size_t i0) const
   {

#ifdef _STATIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, dim0);
#endif

      return this->_data[i0];
   }

//--Assign to the entire array

   template<typename T2>
   Static& operator=(const T2& val) {
      this->assign_constant(val);
      return *this;
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:
   static const size_t _size = dim0;    // Total size of the array

};


/*******************************************************************************
 *
 * External functions
 *
 ******************************************************************************/

/*==============================================================================
 *
 * Routine: operator << for Static class
 *
 * Purpose
 * =======
 *
 *   Writes out a static array.
 *
 * I/O
 * ===
 *
 *   os                 - (I) output stream to write to
 *   A                  - (I) base of static array to write
 *   return             - (O) updated stream
 *
 *============================================================================*/

template <typename T, size_t Size>
std::ostream &operator<<(std::ostream& os, const Static_Base<T, Size> &A)

{

   const T *p = A.begin();
   const T *const pEnd = A.end();
   if ( p != pEnd ) os << *p;
   for ( ++p; p != pEnd; ++p ) os << ' ' << *p;
   return os;

}

}  // End of namespace Array

#endif  // _STATIC_ARRAY_INCLUDED
