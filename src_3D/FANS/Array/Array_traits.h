#ifndef _ARRAY_TRAITS_INCLUDED
#define _ARRAY_TRAITS_INCLUDED


/*******************************************************************************
 *
 * Traits for array classes
 *
 ******************************************************************************/


namespace Array
{


/*==============================================================================
 * Enumerations related to traits
 *============================================================================*/

//--Array size

enum ArraySizeTr_t {
   Small,
   Large
};


//--Data types

enum DataTypeTr_t {
   Primitive,
   Complex
};


/*==============================================================================
 *
 * struct: Traits
 *
 * Purpose
 * =======
 * 
 *   Defines array traits based on dim
 *
 *   size              -- [ Small | Large ] This defines how the arrays are
 *                        traversed.  Small arrays are fully unrolled.  Large
 *                        arrays use a standard loop.
 *
 * Notes
 * =====
 *
 *   - Small arrays have dim 1-60.  The alpha compiler cxx will abort above 65
 *     recursive levels.
 *   - The alpha compiler cxx will not obtain much benefit from unrolling above
 *     50 levels.  Itaniums benefit above 100.  Athlons using gcc benefit up to
 *     about 90.  (Results - 03 Nov 2005)
 *
 *============================================================================*/

template<size_t dim> struct Traits {
   static const ArraySizeTr_t size = Large;
};
template<> struct Traits< 1> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 2> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 3> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 4> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 5> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 6> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 7> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 8> { static const ArraySizeTr_t size = Small; };
template<> struct Traits< 9> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<10> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<11> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<12> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<13> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<14> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<15> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<16> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<17> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<18> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<19> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<20> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<21> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<22> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<23> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<24> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<25> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<26> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<27> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<28> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<29> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<30> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<31> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<32> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<33> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<34> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<35> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<36> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<37> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<38> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<39> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<40> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<41> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<42> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<43> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<44> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<45> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<46> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<47> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<48> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<49> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<50> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<51> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<52> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<53> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<54> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<55> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<56> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<57> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<58> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<59> { static const ArraySizeTr_t size = Small; };
template<> struct Traits<60> { static const ArraySizeTr_t size = Small; };


/*==============================================================================
 *
 * struct: TypeTr
 *
 * Purpose
 * =======
 *
 *   Provides information about a type
 *
 *   isClass           -- [ T | F ] test whether a type is a class or not
 *
 * Notes
 * =====
 *
 *   - The class test used the SFINAE principle as described in 15.2.2
 *     "Determining Class Types" in Vandevoorde and Josuttis "C++ Templates"
 *     book.
 *
 *============================================================================*/

template <typename T>
class TypeTr
{
  private:
   typedef char One;
   typedef struct { char a[2]; } Two;
   template <typename C> static One test(int C::*);
   template <typename C> static Two test(...);
  public:
   enum { IsClass = sizeof(TypeTr<T>::template test<T>(0)) == 1 };
};

}  // End of namespace array

#endif  // _ARRAY_TRAITS_INCLUDED
