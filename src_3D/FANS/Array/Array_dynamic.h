#ifndef _DYNAMIC_ARRAY_INCLUDED
#define _DYNAMIC_ARRAY_INCLUDED


/*==============================================================================
 * Included Files
 *============================================================================*/

//----- Standard Library -----//

#include<functional>
#include<iostream>
#include<new>
#include<stdexcept>
#include<string>


namespace Array
{


/*==============================================================================
 * Forward declarations of sister array types
 *============================================================================*/

//--Static

template <typename T, size_t dim0, size_t dim1, size_t dim2, size_t dim3,
   size_t dim4, size_t dim5, size_t dim6>
class Static;
template <typename T, size_t Size>
class Static_Base;

//--Reference

template <typename T, int Rank>
class Reference;
template <typename T>
class Reference_Base;


/*******************************************************************************
 *
 * class Dynamic
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array (an array for which the dimensions are not known
 *   at compiler time) with up to 7 dimensions by indexing a 1 dimensional
 *   array.  Data space for the arrays may be obtained from "new" or one can
 *    provide a raw address to allocate the data at.
 *
 * Constructors
 * ============
 *
 *   The class takes two template arguments.  This first is the type of data
 *   in the array and the second is the rank of the array.
 *
 *   Dynamic(const size_t dim0 ... dim6)
 *                     -- creates array using new
 *     dim0 - dim6        (I) dimensions of the array (# dimensions = rank)
 *
 *   Dynamic(const size_t dim0 ... dim6, void *const addr)
 *                     -- creates array at address given by addr.
 *     dim0 - dim6        (I) dimensions of the array (# dimensions = rank)
 *     addr               (I) address to create array
 *
 *   Dynamic(const Dynamic &A)
 *                     -- copy constructor using new
 *     A                  (I) array to copy
 *
 *   Dynamic(const Dynamic &A, void *const addr)
 *                     -- copy constructor the creates array at address given
 *                        by addr.
 *     A                  (I) array to copy
 *     addr               (I) address to create array
 *
 *   Dynamic& operator=(const Dynamic &A)
 *                     -- assignment operator.  The two arrays have the same
 *                        type and rank but possibly different dimensions.  Data
 *                        will be copied until one of the two runs out of size.
 *     A                  (I) array to copy
 *     return             (O) updated reference
 *
 *   Dynamic& operator=(const Dynamic<T2, Rank> &A)
 *                     -- assignment operator from a different dynamic array.
 *                        Data will be copied until one of the two runs out of
 *                        size.
 *     A                  (I) array to copy
 *     return             (O) updated reference
 *
 *   Dynamic& operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4,
 *                      SDim5, SDim6> &A)
 *                     -- assignment operator from a static array.  Data will be
 *                        copied until one of the two runs out of size.
 *     A                  (I) array to copy
 *     return             (O) updated reference
 *
 * Member Functions
 * ================
 *
 *   void allocate(const size_t dim0 ... dim6)
 *                     -- allocates the array using new
 *     dim0 - dim6        (I) dimensions of the array (# dimensions = rank)
 *
 *   void allocate(const size_t dim0 ... dim6, void *const addr)
 *                     -- creates array at address given by addr.
 *     dim0 - dim6        (I) dimensions of the array (# dimensions = rank)
 *     addr               (I) address to create array
 *
 *   void deallocate() -- releases memory allocated via new
 *
 *   void reference(Dynamic &A)
 *                     -- Reference the data of array A.  Array A must have the
 *                        same rank as this array.  Some notes:
 *                        - If this array is allocated, it will be deallocated.
 *                        - The dimensions of this array is set to that of A.
 *                        - Deallocation of this array will not affect the data.
 *                        - Deallocate will reset the size of this array so it
 *                          can no longer be used (enable debugging).
 *                        - Do not deallocate A and then index this array!
 *
 *   void reshape(const size_t dim0 ... dim6, T *const addr)
 *                     -- Reference the data at addr and reshape to the rank
 *                        of this array.  Some notes:
 *                        - This is the same as allocating on raw memory except
 *                          constructors and destructors are never called for
 *                          the elements.  It is expected that this array will
 *                          reference the data of another array (which was
 *                          properly allocated).
 *                        - If this array is allocated, it will be deallocated.
 *                        - Deallocation of this array will not affect the data.
 *                        - Deallocate will reset the size of this array so it
 *                          can no longer be used (enable debugging).
 *                        - Do not deallocate the data and then index this
 *                          array!
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
 *   boot is_allocated()
 *                     -- Return boolean indicating if array is allocated
 *     return             (O) T - allocated, F - not allocated
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
 *   operator=(const T2 &val)
 *                     -- assigns "val" to each index in the array (note that
 *                        "val" may be of a different type)
 *     val                (I) value to assign
 *
 *   ostream& operator<<(ostream& os, const Dynamic_Base<T>& A)
 *                     -- writes out the array elements in storage order with
 *                        single space between each element (assumes << defined
 *                        for T)#endif//ARRAY_INCLUDED

 *     os                 (I) output stream
 *     A                  (I) base part of dynamic array
 *     return             (O) output stream
 *
 * Data Members
 * ============
 *
 *   size              -- Overall size of the array
 *
 * Notes
 * =====
 *
 *   - Uses C++ storage order
 *   - Define _DYNAMIC_ARRAY_CHECK_BOUNDS to enable array bounds checking.  If
 *     a bound is exceeded, the value of the incorrect index is printed and the
 *     program aborts.
 *   - There is a base class which implements some general functions.  The
 *     general template class defines arrays of unsupported rank (it aborts).
 *     All supported arrays have specialized template classes.
 *   - The Dynamic class can be defined using the default constructor.
 *     If this is done, memory must be allocated using an allocate function
 *     before the array is used.  There is NO checking for this!
 *
 * Exceptions
 * ==========
 *
 *   std::range_error
 *   std::runtime_error
 *   std::out_of_range
 *
 * Programmer
 * ==========
 *
 *   S. Guzik         15 October 2004     Original code
 *   S. Guzik         14 June 2005        - reorganization such that the base
 *                                        class contains more work.
 *                                        - allocate functions can now be used
 *                                        if actual size of the array is not
 *                                        known during definition.
 *   S. Guzik         08 September 2005   - modified release memory so that raw
 *                                        memory is not deleted, only the
 *                                        elements are deconstructed.
 *
 * Example
 * =======
 *
 *------------------------------------------------------------------------------

#include "Array_dynamic.h"

using std::cout;
using std::endl;

int main() {

//--Common declarations

   int ix1 = 3;
   int ix2 = 4;

//--Create a 2D dynamic array using new

   Array::Dynamic<double, 2> a(ix1, ix2);
   a = 2.;
   a(2, 2) = 4.;
   cout << "sum: " << a.sum() << endl;
   cout << "max: " << a.max() << endl;
   cout << "min: " << a.min() << endl;
   Array::Dynamic<double, 2> b(a);
   cout << "b[2][1]: " << b(2, 1) << " b[2][2]: " << b(2, 2) << endl;
   cout << "All of b: " << b << endl;

}

 *------------------------------------------------------------------------------
 *
 ******************************************************************************/


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic_Base
 *
 * Purpose
 * =======
 *
 *   Provides some base data and functions used by all of the Dynamic
 *   classes.
 *
 ******************************************************************************/

template<typename T>
class Dynamic_Base

{


/*==============================================================================
 * Public types
 *============================================================================*/

  public:

   typedef T value_type;

   enum alloc_type_t {
      alloc_new,
      alloc_raw
   };


/*==============================================================================
 * Public member functions
 *============================================================================*/

//--Max

   const T &max() const
   { 
      return Math::Array1D::maxval(_data, _size, std::less<T>());
   }

   // With user defined less than predicate
   template<typename Cmp>
   const T &max(const Cmp &cmp) const
   {
      return Math::Array1D::maxval(_data, _size, cmp);
   }

//--Min

   const T &min() const
   { 
      return Math::Array1D::minval(_data, _size, std::less<T>());
   }

   // With user defined less than predicate
   template<typename Cmp>

   const T &min(const Cmp &cmp) const
   {
      return Math::Array1D::minval(_data, _size, cmp);
   }

//--Sum

   T sum() const { return Math::Array1D::sum(_data, _size); }

//--Product

   T product() const { return Math::Array1D::product(_data, _size); }

//--Accumulate using user defined accumulation functor

   template<typename Accum>
   T accumulate(const Accum &accum) const
   {
      return Math::Array1D::accumulate(_data, _size, accum);
   }

//--Return pointer to beginning of array

   T *begin() { return _data; }
   const T *begin() const { return _data; }

//--Return pointer to end of array

   T *end() { return _data + _size; }
   const T *end() const { return _data + _size; }

//--Return size of the array

   size_t size() const { return _size; }
   
//--Return allocation status

   bool is_allocated() const { return ( _data ) ? true : false; }

//--Deallocate memory

   void deallocate() {
      release_memory();
      _data = 0;
      _size = 0;
   }


/*==============================================================================
 * Protected constructors and destructors
 *============================================================================*/

  protected:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  Allocation routines must be called before the
//--Dynamic is used.

   Dynamic_Base()
      :
      _data(0), _size(0)
   { }

//--Constructor using new

   Dynamic_Base(const size_t sizeIn)
      :
      _data(0), _size(sizeIn)
   {
      allocate_base();
   }

//--Construct on raw address

   Dynamic_Base(const size_t sizeIn, void *const addr)
      :
      _data(0), _size(sizeIn)
   {
      allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Copy constructor using new

   Dynamic_Base(const Dynamic_Base &A)
      :
      _data(0), _size(A._size)
   {
      if ( A.is_allocated() ) {
         allocate_base();
         assign_base(A);
      }
   }

//--Copy constructor using raw address

   Dynamic_Base(const Dynamic_Base &A, void *const addr)
      :
      _data(0), _size(A._size)
   {
      if ( A.is_allocated() ) {
         allocate_base(addr);
         assign_base(A);
      }
   }

/*--------------------------------------------------------------------*
 * Destructor
 *--------------------------------------------------------------------*/

   ~Dynamic_Base()
   {
      release_memory();
   }


/*==============================================================================
 * Protected member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assignment of base
 *--------------------------------------------------------------------*/

   // May be same or different types
   template<typename T2>
   void assign_base(const Dynamic_Base<T2> &A)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      check_allocated();
#endif

      T* p = _data;
      const T2* a = A.begin();
      for ( size_t n = std::min(_size, A._size); n--; ) *p++ = *a++;
   }

   // From a static array
   template <typename T2, size_t StaticSize>
   void assign_base(const Static_Base<T2, StaticSize> &A)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      check_allocated();
#endif

      T *p = _data;
      const T2 *a = A.begin();
      for ( size_t n = std::min(_size, StaticSize); n--; ) *p++ = *a++;
   }

   // From a reference array
   template <typename T2>
   void assign_base(const Reference_Base<T2> &A)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      check_allocated();
#endif

      T *p = _data;
      const T2 *a = A.begin();
      for ( size_t n = std::min(_size, A.size()) ; n-- ; ) *p++ = *a++;
   }

/*--------------------------------------------------------------------*
 * Assignment of constant to base
 *--------------------------------------------------------------------*/

   template<typename T2>
   void assign_constant(const T2 &val)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      check_allocated();
#endif

      T* p = _data;
      for ( size_t n = _size ; n-- ; ) *p++ = val;
   }

/*--------------------------------------------------------------------*
 * Allocate the base memory
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate_base()
   {
      if ( _data ) release_memory();
      _data = new(std::nothrow) T[_size];
      if ( _data == 0 ) {
         throw std::runtime_error("Allocation of memory failed for class "
                                  "Array::Dynamic");
      }
      allocType = alloc_new;
   }

//--Allocate on raw address

   void allocate_base(void *const addr)
   {
      if ( _data ) release_memory();
      _data = AllocRawPolicy<T>::eval(_size, addr);
      allocType = alloc_raw;
   }


/*==============================================================================
 * Private constructors
 *============================================================================*/

  private:

//--Standard assignmnet not permitted

   Dynamic_Base& operator=(const Dynamic_Base&);


/*==============================================================================
 * Private member functions
 *============================================================================*/

   void release_memory() {
      switch ( allocType ) {
      case alloc_new:
         delete[] _data;
         break;
      case alloc_raw:
         // Deconstruct the elements if they have class type
         ReleaseRawPolicy<T>::eval(_size, _data);
         break;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  protected:

   T* _data;                            // Data for the array
   size_t _size;                        // Overall size of the array

  private:

   alloc_type_t allocType;              // Type of allocation


/*==============================================================================
 * Tools for debugging
 *============================================================================*/

  protected:

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
   void check_allocated() const
   {
      if ( _data == 0 ) {
         throw std::domain_error("Array not allocated\n"
                                 "Error in class Array::Reference");
      }
   }
#endif

#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
   void check_bounds(const int dim, const size_t index, const size_t iMax) const
   {
      if ( index < 0 || index >= iMax ) {
         char oText[128];
         std::sprintf(oText, "Index %d is outside array bounds:\n"
                      "Index value: %d\n"
                      "Valid range: 0 : %d\n"
                      "Error in class Array::Dynamic", dim, index, iMax-1);
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
      return MPIuT.Create_contiguous(_size);
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


/*==============================================================================
 * Friends
 *============================================================================*/

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - unsupported rank
 *
 * Purpose
 * =======
 *
 *   General Dynamic template class.  This class aborts because any array
 *   with a rank that has not been specialized is unsupported.
 *
 ******************************************************************************/


template<typename T, int Rank>
class Dynamic
   
{

  public:

   Dynamic() {
      throw std::range_error("Arrays greater than rank 7 not supported for "
                             "class Array::Dynamic");
   }

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 7
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 7.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 7>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4, const size_t dim5,
           const size_t dim6)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4*dim5*dim6),
      stride5(dim6),
      stride4(stride5*dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4, const size_t dim5,
           const size_t dim6, void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4*dim5*dim6, addr),
      stride5(dim6),
      stride4(stride5*dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride5(A.stride5),
      stride4(A.stride4),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride5(A.stride5),
      stride4(A.stride4),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic& operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4, const size_t i5,
                 const size_t i6)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3/stride4);
      this->check_bounds(5, i5, stride4/stride5);
      this->check_bounds(6, i6, stride5);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5*stride5 + i6];
   }

   const T &operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4, const size_t i5,
                       const size_t i6) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3/stride4);
      this->check_bounds(5, i5, stride4/stride5);
      this->check_bounds(6, i6, stride5);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5*stride5 + i6];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4, const size_t dim5,
                 const size_t dim6)
   {
      stride5 = dim6;
      stride4 = stride5*dim5;
      stride3 = stride4*dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4, const size_t dim5,
                 const size_t dim6, void *const addr)
   {
      stride5 = dim6;
      stride4 = stride5*dim5;
      stride3 = stride4*dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0/stride1;
         break;
      case 3:
         return stride1/stride2;
         break;
      case 4:
         return stride2/stride3;
         break;
      case 5:
         return stride3/stride4;
         break;
      case 6:
         return stride4/stride5;
         break;
      case 7:
         return stride5;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride5;                      // Stride for index 5
   size_t stride4;                      // Stride for index 4
   size_t stride3;                      // Stride for index 3
   size_t stride2;                      // Stride for index 2
   size_t stride1;                      // Stride for index 1
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 7>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 6
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 6.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 6>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4, const size_t dim5)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4*dim5),
      stride4(dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4, const size_t dim5,
           void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4*dim5, addr),
      stride4(dim5),
      stride3(stride4*dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride4(A.stride4),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride4(A.stride4),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A) {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank> &A) {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4, const size_t i5)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3/stride4);
      this->check_bounds(5, i5, stride4);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5];
   }

   const T &operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4, const size_t i5) const 
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3/stride4);
      this->check_bounds(5, i5, stride4);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4*stride4 + i5];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4, const size_t dim5)
   {
      stride4 = dim5;
      stride3 = stride4*dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4, const size_t dim5,
                 void *const addr)
   {
      stride4 = dim5;
      stride3 = stride4*dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0/stride1;
         break;
      case 3:
         return stride1/stride2;
         break;
      case 4:
         return stride2/stride3;
         break;
      case 5:
         return stride3/stride4;
         break;
      case 6:
         return stride4;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride4;                      // Stride for index 4
   size_t stride3;                      // Stride for index 3
   size_t stride2;                      // Stride for index 2
   size_t stride1;                      // Stride for index 1
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 6>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 5
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 5.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 5>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4),
      stride3(dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, const size_t dim4, void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3*dim4, addr),
      stride3(dim4),
      stride2(stride3*dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride3(A.stride3),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3, const size_t i4)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4];
   }

   const T &operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3, const size_t i4) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2/stride3);
      this->check_bounds(4, i4, stride3);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3*stride3 +
                         i4];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4)
   {
      stride3 = dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, const size_t dim4, void *const addr)
   {
      stride3 = dim4;
      stride2 = stride3*dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0/stride1;
         break;
      case 3:
         return stride1/stride2;
         break;
      case 4:
         return stride2/stride3;
         break;
      case 5:
         return stride3;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride3;                      // Stride for index 3
   size_t stride2;                      // Stride for index 2
   size_t stride1;                      // Stride for index 1
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 5>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 4
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 4.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 4>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3),
      stride2(dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           const size_t dim3, void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1*dim2*dim3, addr),
      stride2(dim3),
      stride1(stride2*dim2),
      stride0(stride1*dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride2(A.stride2),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank>& A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1, const size_t i2,
                 const size_t i3)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3];
   }

   const T &operator()(const size_t i0, const size_t i1, const size_t i2,
                       const size_t i3) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1/stride2);
      this->check_bounds(3, i3, stride2);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2*stride2 + i3];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3)
   {
      stride2 = dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 const size_t dim3, void *const addr)
   {
      stride2 = dim3;
      stride1 = stride2*dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0/stride1;
         break;
      case 3:
         return stride1/stride2;
         break;
      case 4:
         return stride2;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride2;                      // Stride for index 2
   size_t stride1;                      // Stride for index 1
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 4>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 3
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 3.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 3>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2)
      :
      Dynamic_Base<T>(dim0*dim1*dim2),
      stride1(dim2),
      stride0(stride1*dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, const size_t dim2,
           void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1*dim2, addr),
      stride1(dim2),
      stride0(stride1*dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride1(A.stride1),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1, const size_t i2)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2];
   }

   const T &operator()(const size_t i0, const size_t i1, const size_t i2) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0/stride1);
      this->check_bounds(2, i2, stride1);
#endif

      return this->_data[i0*stride0 + i1*stride1 + i2];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2)
   {
      stride1 = dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, const size_t dim2,
                 void *const addr)
   {
      stride1 = dim2;
      stride0 = stride1*dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0/stride1;
         break;
      case 3:
         return stride1;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride1;                      // Stride for index 1
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 3>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 2
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 2.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 2>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic() :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0, const size_t dim1)
      :
      Dynamic_Base<T>(dim0*dim1),
      stride0(dim1)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, const size_t dim1, void *const addr)
      :
      Dynamic_Base<T>(dim0*dim1, addr),
      stride0(dim1)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A),
      stride0(A.stride0)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr),
      stride0(A.stride0)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank>& A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0, const size_t i1)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0);
#endif

      return this->_data[i0*stride0 + i1];
   }

   const T &operator()(const size_t i0, const size_t i1) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size/stride0);
      this->check_bounds(1, i1, stride0);
#endif

      return this->_data[i0*stride0 + i1];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0, const size_t dim1)
   {
      stride0 = dim1;
      this->_size = stride0*dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, const size_t dim1, void *const addr)
   {
      stride0 = dim1;
      this->_size = stride0*dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
         return this->_size;
         break;
      case 1:
         return this->_size/stride0;
         break;
      case 2:
         return stride0;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   size_t stride0;                      // Stride for index 0


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class Reference<T, 2>;

/*    friend LAY::Layout LAY::array_layout( */
/*       const char aType, const unsigned int aDim, const int nElem, */
/*       const int elemSize, int* counts, Layout** layoutTypes, MPI::Aint* displs, */
/*       const char *const arrayName, const char *const dataName); */

};


/*******************************************************************************
 *\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 *******************************************************************************
 *
 * class Dynamic - rank 1
 *
 * Purpose
 * =======
 *
 *   Implements a dynamic array of rank 1.  This is a specialized template class
 *   of the general unsupported rank template class.
 *
 ******************************************************************************/

template<typename T>
class Dynamic<T, 1>: public Dynamic_Base<T>

{


/*==============================================================================
 * Public constructors and destructors
 *============================================================================*/

  public:

/*--------------------------------------------------------------------*
 * Standard constructors
 *--------------------------------------------------------------------*/

//--Default constructor.  The allocate function must be called before the
//--dynamic array is used.

   Dynamic()
      :
      Dynamic_Base<T>()
   { }

//--Construct using new

   Dynamic(const size_t dim0)
      :
      Dynamic_Base<T>(dim0)
   { }

//--Construct on raw address

   Dynamic(const size_t dim0, void *const addr)
      :
      Dynamic_Base<T>(dim0, addr)
   { }

/*--------------------------------------------------------------------*
 * Copy constructors
 *--------------------------------------------------------------------*/

//--Construct using new

   Dynamic(const Dynamic &A)
      :
      Dynamic_Base<T>(A)
   { }

//--Construct on raw address

   Dynamic(const Dynamic &A, void *const addr)
      :
      Dynamic_Base<T>(A, addr)
   { }

/*--------------------------------------------------------------------*
 * Assignment operator
 *--------------------------------------------------------------------*/

//--Use base class assignment operator

   // Same type
   Dynamic &operator=(const Dynamic &A)
   {
      if ( &A != this ) this->assign_base(A);
      return *this;
   }

   // Different type and/or rank
   template<typename T2, int Rank>
   Dynamic &operator=(const Dynamic<T2, Rank>& A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a static array
   template <typename T2, size_t SDim0, size_t SDim1, size_t SDim2,
      size_t SDim3, size_t SDim4, size_t SDim5, size_t SDim6>
   Dynamic &operator=(const Static<T2, SDim0, SDim1, SDim2, SDim3, SDim4, SDim5,
                      SDim6> &A)
   {
      this->assign_base(A);
      return *this;
   }

   // From a reference array
   template<typename T2, int Rank>
   Dynamic &operator=(const Reference<T2, Rank> &A)
   {
      this->assign_base(A);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Use synthesized destructor
 *--------------------------------------------------------------------*/


/*==============================================================================
 * Public operators
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Assign constant to the array
 *--------------------------------------------------------------------*/

   template<typename T2>
   Dynamic &operator=(const T2 &val)
   {
      this->assign_constant(val);
      return *this;
   }

/*--------------------------------------------------------------------*
 * Index the array
 *--------------------------------------------------------------------*/

   T &operator()(const size_t i0)
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size);
#endif

      return this->_data[i0];
   }

   const T &operator()(const size_t i0) const
   {

#ifdef _DYNAMIC_ARRAY_CHECK_ALLOCATED
      this->check_allocated();
#endif
#ifdef _DYNAMIC_ARRAY_CHECK_BOUNDS
      this->check_bounds(0, i0, this->_size);
#endif

      return this->_data[i0];
   }


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Allocate the memory (if not performed during construction)
 *--------------------------------------------------------------------*/

//--Allocate using new

   void allocate(const size_t dim0)
   {
      this->_size = dim0;
      this->allocate_base();
   }

//--Allocate on raw address

   void allocate(const size_t dim0, void *const addr)
   {
      this->_size = dim0;
      this->allocate_base(addr);
   }

/*--------------------------------------------------------------------*
 * Return size of dimensions
 *--------------------------------------------------------------------*/

   size_t dim_size(const size_t dim = 0) const
   {
      switch ( dim ) {
      case 0:
      case 1:
         return this->_size;
         break;
      default:
         return 0;
      }
   }


/*==============================================================================
 * Data members
 *============================================================================*/

  private:

};


/*******************************************************************************
 *
 * External functions
 *
 ******************************************************************************/

/*==============================================================================
 *
 * Routine: operator << for Dynamic_Base class
 *
 * Purpose
 * =======
 *
 *   Writes out a dynamic array.
 *
 * I/O
 * ===
 *
 *   os                 - (I) output stream to write to
 *   Dynamic            - (I) dynamic array to write
 *   return             - (O) updated stream
 *
 *============================================================================*/

template<typename T>
std::ostream &operator<<(std::ostream &os, const Dynamic_Base<T> &A)

{

   const T *p = A.begin();
   const T *const pEnd = A.end();
   if ( p != pEnd ) os << *p;
   for ( ++p; p != pEnd; ++p ) os << ' ' << *p;
   return os;

}

}  // End of namespace Array

#endif  // _DYNAMIC_ARRAY_INCLUDED
