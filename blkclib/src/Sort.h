#ifndef _SORT_INCLUDED
#define _SORT_INCLUDED

/*******************************************************************************

   Copyright (C) 2007 [--DRAFT VERSION--NOT FOR DISTRIBUTION--]

   This file is a part of the Block Connectivity library 'libblkc'

   'libblkc' is free software; you can redistribute it and/or modify it under
   the terms of the GNU General Public License as published by the Free Software
   Foundation; either version 3 of the License, or (at your option) any later
   version.

   'libblkc' is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
   more details.

   You should have received a copy of the GNU General Public License along with
   this program.  If not, see <http://www.gnu.org/licenses/>.

*******************************************************************************/


/*==============================================================================
 * Included Files
 *============================================================================*/

//----- Standard Library -----//

#include<functional>
#include<limits>
#include<stdexcept>
#include<string>


/*==============================================================================
 * Set file parameters
 *============================================================================*/

const std::string fn_Sort_h("Sort.h");  // Name of this file


namespace Sort
{

const int defaultThreshold = 27;        // Default threshold.  Typically in
                                        // the range of 15-35.


/*==============================================================================
 *
 * Routine: log2
 *
 * Purpose
 * =======
 *
 *   Simple and fast routine for performing log2 of an integer.  The result
 *   is the floor of the log2 of the real.
 *
 * Notes
 * =====
 *
 *   - val must be > 0
 *
 *============================================================================*/

inline unsigned long log2(unsigned long val) 

{

   unsigned long result = std::numeric_limits<unsigned long>::max();
   while ( val ) {
      ++result;
      val >>= 1;
   }
   return result;

}


/*==============================================================================
 * Forward declarations of sort routines
 *============================================================================*/

template<typename T, typename Cmp>
inline void quick(int n, T *array, const int threshold, int depth,
                  const Cmp& cmp); 

template<typename T, typename Cmp>
inline void heap(const int n, T *array, const Cmp& cmp);

template<typename T, typename Cmp>
inline void insertion(const int n, T *array, Cmp& cmp); 


/*==============================================================================
 *
 * Routine: intro
 *
 * Purpose
 * =======
 *
 *   Performs an introspective sort on an array.  Routine quick does all the
 *   work.  This provides setup and it is overloaded.  The first uses standard
 *   < to compare the objects T (as provided by the standard library predicate
 *   std::less<T>).  In the second, the user can provide their own predicate.
 *
 * I/O
 * ===
 *
 *   n                  - (I) number of elements in the array
 *   array              - (I) array to be sorted.  Elements numbered 0 to n-1
 *                            and are contiguous in memory
 *                      - (O) the sorted array
 *   cmp                - (I) binary predicate to test if arg1 < arg2.
 *   threshold          - (I) size of partitions before switching to an
 *                            insertion sort (default 27)
 *
 *============================================================================*/


/*------------------------------------------------------------------------------
 * Version with default predicate
 *----------------------------------------------------------------------------*/

template<typename T>
inline void intro(int n, T* array, const int threshold = defaultThreshold)

{

   if ( threshold < 2 ) throw std::invalid_argument(
      "Threshold must be >= 2 in function intro in file " + fn_Sort_h);

   // The normal depth for a quicksort without a threshold is 2*log2(n) [See
   // Musser, "Instrospective Sorting and Selection Algorithms", RPI.].  For
   // perfect partitioning, this translates to 2*log2(n/t) with threshold t.
   // For some reason, this seems a bit wishful.  Halfway in between is t/2.  So
   // instead we use 2*log2(2*n/t).

   const int depth = 2*log2(std::max(1, 2*int(double(n)/double(threshold))));
   quick(n, array, threshold, depth, std::less<T>());

   // Uncomment this line and comment out the one at the end of routine quick
   // to insertion sort the full array at the end.  Insertion sorting the parts
   // seems to be a tiny bit faster (AMD_64 3500)

//   insertion(n, array);
   
}


/*------------------------------------------------------------------------------
 * Version with user-defined predicate
 *----------------------------------------------------------------------------*/

template<typename T, typename Cmp>
inline void intro(int n, T* array, Cmp cmp,
                  const int threshold = defaultThreshold)

{

   if ( threshold < 2 ) throw std::invalid_argument(
      "Threshold must be >= 2 in function intro in file " + fn_Sort_h);
   const int depth = 2*log2(std::max(1, 2*int(double(n)/double(threshold))));
   quick(n, array, threshold, depth, cmp);
//   insertion(n, array);
   
}


/*==============================================================================
 *
 * Routine: quick
 *
 * Purpose
 * =======
 *
 *   Performs a introspective sort on an array.  This is essentially a
 *   quicksort that resorts to a heapsort when O(N^2) behaviour is detected.
 *
 * I/O
 * ===
 *
 *   n                  - (I) number of elements in the array
 *   array              - (I) array to be sorted.  Elements numbered 0 to n-1
 *                            and are contiguous in memory
 *                      - (O) the sorted array
 *   threshold          - (I) size of partitions before switching to an
 *                            insertion sort
 *   depthLevel         - (I) current level of recursion
 *
 *============================================================================*/


template<typename T, typename Cmp>
inline void quick(int n, T *array, const int threshold, int depth,
                  const Cmp& cmp)

{

   const int l = 0;                     // Left-most point
   const int p = 1;                     // Stores partition value
   int r = n-1;                         // Right-most point

//--Loop and recurse until size of partition is less than threshold

   while ( n > threshold ) {

//--Check for O(N^2) behaviour and use heapsort if detected

      if ( depth == 0 ) {
         heap(n, array, cmp);
         return;
      }
      --depth;

      int k = n/2;                      // Middle point
      T partition;                      // Partitioning element

//--Use median of 3 to find partition.  Arrange l, k, and r elements into l, p,
//--and r with order f(l) < f(p) < f(r).

      if ( cmp(array[l], array[k]) ) {
         if ( cmp(array[k], array[r]) ) {
            // k is middle element, put in spot p
            partition = array[k];
            array[k] = array[p];
         }
         else if ( cmp(array[l], array[r]) ) {
            // r is middle element, put in spot p.  Put k in r and p in k.
            partition = array[r];
            array[r] = array[k];
            array[k] = array[p];
         }
         else {
            // l is middle element, put in spot p.  Put k in r and r in l and
            // p in k.
            partition = array[l];
            array[l] = array[r];
            array[r] = array[k];
            array[k] = array[p];
         }
      }
      else {
         if ( cmp(array[r], array[k]) ) {
            // k is middle element, swap with spot p.  Swap l and r.
            partition = array[l];
            array[l] = array[r];
            array[r] = partition;
            partition = array[k];
            array[k] = array[p];
         }
         else if ( cmp(array[l], array[r]) ) {
            // l is middle element, put in spot p.  Put k in l and p in k.
            partition = array[l];
            array[l] = array[k];
            array[k] = array[p];
         }
         else {
            // r is middle element, put in spot p.  Put k in l and l in r and
            // p in k.
            partition = array[r];
            array[r] = array[l];
            array[l] = array[k];
            array[k] = array[p];
         }
      }
      array[p] = partition;  // Required stopper for backwards pointer (even
                             // though this is always retained in "partition").

//--Elements l, p, and r are now in increasing order.  p is the partition.  k is
//--no longer used.

      T *forward = &array[p];   // Start at p
      T *backward = &array[r];  // Start at r

//--Partition around element p

      while ( true ) {
         do ++forward; while ( cmp(*forward, partition) );
         do --backward; while ( cmp(partition, *backward) );

         if ( forward > backward ) break;  // Pointers crossed

         // Swap elements pointed to by forward and backward

         T temp = *forward;
         *forward = *backward;
         *backward = temp;
      }

      // Insert the partitioning element into backward
      array[p] = *backward;
      *backward = partition;

//--Recurse on upper side of partition and iterate on lower side

      r = backward - array - 1;
      const int lu = r + 2;
      quick(n - lu, &array[lu], threshold, depth, cmp);
      n = r+1;

   }

//--Finish the partition with an insertion sort (comment out if insertion to
//--be performed on entire array at the end).

   insertion(n, array, cmp);

}


/*==============================================================================
 *
 * Routine: heap
 *
 * Purpose
 * =======
 *
 *   Performs a heap sort on an array
 *
 * I/O
 * ===
 *
 *   n                  - (I) Size of the array
 *   array              - (I) first index of unsorted array
 *                      - (O) sorted array
 *
 * Notes
 * =====
 *
 *   Requires n > 1 (no checking performed)
 *
 *============================================================================*/

template<typename T, typename Cmp>
inline void heap(const int n, T *array, const Cmp& cmp)

{

   T temp;

//--Load the heap

   int r = n-1;
   for ( int l = n/2 - 1 ; l >= 0 ; --l ) {
      temp = array[l];  // Value starting in new position
      int i = l;        // Current position available
      int j = i+i+1;    // First of competitors for new position
      while ( j < n ) {  // The shift down loop
         if ( j < r ) {
            if ( cmp(array[j], array[j+1]) ) ++j;
         }
         if ( cmp(temp, array[j]) ) {  // Move temp down
            array[i] = array[j];
            i = j;
            j += j+1;
         }
         else {  // This is where temp goes so exit loop
            j = n;
         }
      }
      array[i] = temp;
   }

//--Unload the heap

   while ( r != 1 ) {
      temp = array[r];
      array[r] = array[0];
      --r;
      int i = 0;  // Current position available
      int j = 1;  // First of competitors for new position
      while ( j <= r ) {  // The shift down loop
         if ( j < r ) {
            if ( cmp(array[j], array[j+1]) ) ++j;
         }
         if ( cmp(temp, array[j]) ) {  // Move temp down
            array[i] = array[j];
            i = j;
            j += j+1;
         }
         else {  // This is where temp goes so exit loop
            j = n;
         }
      }
      array[i] = temp;
   }
   temp = array[1];
   array[1] = array[0];
   array[0] = temp;

}


/*==============================================================================
 *
 * Routine: insertion
 *
 * Purpose
 * =======
 *
 *   Performs an insertion sort on an array
 *
 * I/O
 * ===
 *
 *   n                  - (I) Size of the array
 *   array              - (I) first index of unsorted array
 *                      - (O) sorted array
 *
 *============================================================================*/

template<typename T, typename Cmp>
inline void insertion(const int n, T *array, Cmp& cmp)

{

   T temp;
   for ( int j = 1 ; j < n ; ++j ) {
      T temp(array[j]);
      int i(j-1);
      while ( cmp(temp, array[i]) ) {
         array[i+1] = array[i];
         --i;
         if ( i < 0 ) break;
      }
      array[i+1] = temp;
   }

}


/*==============================================================================
 *
 * Routine find
 *
 * Purpose
 * =======
 *
 *   Finds index of a given object in a SORTED array of size n.  The object is
 *   assumed to exist and be unique.  If not unique, any valid index may be
 *   returned.  If it doesn't exist, the index of the next lowest value will
 *   probably be found.
 *
 * Notes
 * =====
 *
 *   - The threshold for switching to a sequential search is the same as that
 *   used while sorting to switch to an insertion sort.  The exact setting
 *   probably doesn't matter much.
 *
 *============================================================================*/

template<typename T, typename Cmp>
inline int find(const T &val, const int n, const T *const A, Cmp cmp)
{

   // Binary tree search
   int l = 0;
   int r = n - 1;
   while ( (r - l) > defaultThreshold ) {
      const int m = (l+r)/2;
      const T mVal = A[m];
      if ( cmp(val, mVal) ) r = m;
      else if ( cmp(mVal, val) ) l = m;
      else return m;
   }

   // Sequential search
   for ( int i = l+1 ; i <= r ; ++i ) {
      if ( cmp(val, A[i]) ) return i-1;
   }
   return r;

}

template<typename T>
inline int find(const T &val, const int n, const T *const A) {
   return find(val, n, A, std::less<T>());
}

}  // End of namespace Sort

#endif  // _SORT_INCLUDED
