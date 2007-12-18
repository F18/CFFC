#ifndef _ARRAY_POLICY_INCLUDED
#define _ARRAY_POLICY_INCLUDED


/*******************************************************************************
 *
 * Policies for the array class
 *
 ******************************************************************************/


/*==============================================================================
 * Included Files
 *============================================================================*/

//----- Standard Library -----//

#include <algorithm>
#include <functional>

//----- Template Sources -----//

#include "Array_traits.h"
#include "Math_Array1D.h"


namespace Array
{


/*==============================================================================
 * Policies for dynamic arrays
 *============================================================================*/


/*--------------------------------------------------------------------*
 * Raw address allocation policy
 *--------------------------------------------------------------------*/

   template <typename T, bool IsClass = TypeTr<T>::IsClass>
   struct AllocRawPolicy;

//--If elements are class type, allocate using placement new

   template <typename T>
   struct AllocRawPolicy<T, true>
   {
      static T *eval(const size_t size, void *const addr)
      {
         return new(addr) T[size];
      }
   };

//--If element are not class type, just assign the address

   template <typename T>
   struct AllocRawPolicy<T, false>
   {
      static T *eval(const size_t size, void *const addr)
      {
         return static_cast<T*>(addr);
      }
   };


/*--------------------------------------------------------------------*
 * Raw address release policy
 *--------------------------------------------------------------------*/

   template <typename T, bool IsClass = TypeTr<T>::IsClass>
   struct ReleaseRawPolicy;

//--If elements are class type, invoke destructor on each element

   template <typename T>
   struct ReleaseRawPolicy<T, true>
   {
      static void eval(const size_t size, T *p)
      {
         for ( size_t n = size; n--; ) p++->~T();
      }
   };

//--If elements are not class type, do nothing

   template <typename T>
   struct ReleaseRawPolicy<T, false>
   {
      static void eval(const size_t size, T *p) { }
   };


/*==============================================================================
 * Policies for static arrays
 *============================================================================*/


/*--------------------------------------------------------------------*
 * Assignment policy (array to array)
 *--------------------------------------------------------------------*/

   template<typename T1, typename T2, size_t Dim, ArraySizeTr_t ArraySizeTr>
   struct AssignPolicy;

//--Small arrays

   template<typename T1, typename T2, size_t Dim>
   struct AssignPolicy<T1, T2, Dim, Small> {
      static void eval(T1* p, const T2* a) {
         *p = *a;
         AssignPolicy<T1, T2, Dim-1, Small>::eval(p+1, a+1);
      }
   };

   template<typename T1, typename T2>
   struct AssignPolicy<T1, T2, 1, Small> {
      static void eval(T1* p, const T2* a) { *p = *a; }
   };

//--Large arrays

   template<typename T1, typename T2, size_t Dim>
   struct AssignPolicy<T1, T2, Dim, Large> {
      static void eval(T1* p, const T2* a) {
         for ( size_t n = Dim ; n-- ; ) *p++ = *a++;
      }
   };


/*--------------------------------------------------------------------*
 * Assignment policy (const to array)
 *--------------------------------------------------------------------*/

   template<typename T1, typename T2, size_t Dim, ArraySizeTr_t ArraySizeTr>
   struct AssignConstantPolicy;

//--Small arrays

   template<typename T1, typename T2, size_t Dim>
   struct AssignConstantPolicy<T1, T2, Dim, Small> {
      static void eval(T1* p, const T2& val) {
         *p = val;
         AssignConstantPolicy<T1, T2, Dim-1, Small>::eval(p+1, val);
      }
   };

   template<typename T1, typename T2>
   struct AssignConstantPolicy<T1, T2, 1, Small> {
      static void eval(T1* p, const T2& val) { *p = val; }
   };

//--Large arrays

   template<typename T1, typename T2, size_t Dim>
   struct AssignConstantPolicy<T1, T2, Dim, Large> {
      static void eval(T1* p, const T2& val) {
         for ( size_t n = Dim ; n-- ; ) *p++ = val;
      }
   };


/*--------------------------------------------------------------------*
 * Minimum value
 *--------------------------------------------------------------------*/

   template<typename T, size_t Dim, ArraySizeTr_t ArraySizeTr, typename Cmp>
   struct MinValPolicy;

//--Small arrays

   template<typename T, size_t Dim, typename Cmp>
   struct MinValPolicy<T, Dim, Small, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return std::min(
            *p, MinValPolicy<T, Dim-1, Small, Cmp>::eval(p+1, cmp), cmp);
      }
   };

   template<typename T, typename Cmp>
   struct MinValPolicy<T, 1, Small, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return *p;
      }
   };

//--Large arrays

   template<typename T, size_t Dim, typename Cmp>
   struct MinValPolicy<T, Dim, Large, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return Math::Array1D::minval(p, Dim, cmp);
      }
   };


/*--------------------------------------------------------------------*
 * Maximum value
 *--------------------------------------------------------------------*/

   template<typename T, size_t Dim, ArraySizeTr_t ArraySizeTr, typename Cmp>
   struct MaxValPolicy;

//--Small arrays

   template<typename T, size_t Dim, typename Cmp>
   struct MaxValPolicy<T, Dim, Small, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return std::max(
            *p, MaxValPolicy<T, Dim-1, Small, Cmp>::eval(p+1, cmp), cmp);
      }
   };

   template<typename T, typename Cmp>
   struct MaxValPolicy<T, 1, Small, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return *p;
      }
   };

//--Large arrays

   template<typename T, size_t Dim, typename Cmp>
   struct MaxValPolicy<T, Dim, Large, Cmp> {
      static const T& eval(const T* p, const Cmp& cmp) {
         return Math::Array1D::maxval(p, Dim, cmp);
      }
   };


/*--------------------------------------------------------------------*
 * Accumulation (sum, product, etc.)
 *--------------------------------------------------------------------*/

   template<typename T, size_t Dim, ArraySizeTr_t ArraySizeTr, typename Accum>
   struct AccumulatePolicy;

//--Small arrays

   template<typename T, size_t Dim, typename Accum>
   struct AccumulatePolicy<T, Dim, Small, Accum> {
      static T eval(const T* p, const Accum& accum) {
         return accum(
            *p, AccumulatePolicy<T, Dim-1, Small, Accum>::eval(p+1, accum));
      }
   };

   template<typename T, typename Accum>
   struct AccumulatePolicy<T, 1, Small, Accum> {
      static const T& eval(const T* p, const Accum& accum) {
         return *p;
      }
   };

//--Large arrays

   template<typename T, size_t Dim, typename Accum>
   struct AccumulatePolicy<T, Dim, Large, Accum> {
      static T eval(const T* p, const Accum& accum) {
         return Math::Array1D::accumulate(p, Dim, accum);
      }
   };

}  // End of namespace Array

#endif  // _ARRAY_POLICY_INCLUDED
