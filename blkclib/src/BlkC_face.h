#ifndef _BLKC_FACE_INCLUDED
#define _BLKC_FACE_INCLUDED

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

--------------------------------------------------------------------------------

   Adapted from GMSH source code - Copyright (C) 1997-2007 C. Geuzaine,
   J.-F. Remacle and licensed under GNU General Public License version 2.

*******************************************************************************/

#include <functional>

#include "hash.h"
#include "BlkC_edge.h"

namespace BlkC
{


/*******************************************************************************
 *
 * Class Face
 *
 * Purpose
 * =======
 *
 *   A comparable representation of a block face.
 *
 ******************************************************************************/

class Face
{
  private:
   int _v[4];
   int _si[4];                          // sorted indices to _v

  public:
   Face(const int v0, const int v1, const int v2, const int v3)
   {
      _v[0] = v0; _v[1] = v1; _v[2] = v2; _v[3] = v3;
      sort_vertex();
   }
   Face(const int *const v)
   {
      _v[0] = v[0]; _v[1] = v[1]; _v[2] = v[2]; _v[3] = v[3];
      sort_vertex();
   }
   // Use sythesized copy, assignment, and destructor

   unsigned int get_num_vertices() const { return 4; }
   int get_vertex(const int i) const { return _v[i]; }
   int get_sorted_vertex(const int i) const { return _v[_si[i]]; }
   void get_ordered_vertices(int *const v) const
   {
      v[0] = get_sorted_vertex(0);
      v[1] = get_sorted_vertex(1);
      v[2] = get_sorted_vertex(2);
      v[3] = get_sorted_vertex(3);
   }
   bool has_edge(const Edge &edge) 
   {
      int iE = 0;
      if ( edge.get_sorted_vertex(iE) == get_sorted_vertex(0) ) ++iE;
      if ( edge.get_sorted_vertex(iE) == get_sorted_vertex(1) ) ++iE;
      if ( iE == 2 ) return true;
      if ( edge.get_sorted_vertex(iE) == get_sorted_vertex(2) ) ++iE;
      if ( iE == 2 ) return true;
      if ( edge.get_sorted_vertex(iE) == get_sorted_vertex(3) ) ++iE;
      if ( iE == 2 ) return true;
      return false;
   }
   bool has_vertex(const int vert)
   {
      if ( vert == get_vertex(0) || vert == get_vertex(1) ||
           vert == get_vertex(2) || vert == get_vertex(3) ) return true;
      return false;
   }

  private:
   void sort_vertex()
   {
      // This is simply an unrolled insertion sort (hopefully fast).
      if ( _v[1] < _v[0] ) {
         _si[0] = 1;
         _si[1] = 0;
      }
      else {
         _si[0] = 0;
         _si[1] = 1;
      }
      if ( _v[2] < _v[_si[1]] ) {
         _si[2] = _si[1];
         if ( _v[2] < _v[_si[0]] ) {
            _si[1] = _si[0];
            _si[0] = 2;
         }
         else
            _si[1] = 2;
      }
      else
         _si[2] = 2;
      if ( _v[3] < _v[_si[2]] ) {
         _si[3] = _si[2];
         if ( _v[3] < _v[_si[1]] ) {
            _si[2] = _si[1];
            if ( _v[3] < _v[_si[0]] ) {
               _si[1] = _si[0];
               _si[0] = 3;
            }
            else
               _si[1] = 3;
         }
         else
            _si[2] = 3;
      }
      else
         _si[3] = 3;
   }
};

//--Operators for comparing faces

bool operator==(const Face &f1, const Face &f2) 
{
   return ( f1.get_sorted_vertex(0) == f2.get_sorted_vertex(0) &&
            f1.get_sorted_vertex(1) == f2.get_sorted_vertex(1) &&
            f1.get_sorted_vertex(2) == f2.get_sorted_vertex(2) &&
            f1.get_sorted_vertex(3) == f2.get_sorted_vertex(3) );
}

//--The following function objects compare the addresses of the mesh vertices.
//--Equal, Less, and a Hash are defined.

struct Equal_Face : public std::binary_function<Face, Face, bool> {
   bool operator()(const Face &f1, const Face &f2) const
   {
      return ( f1 == f2 );
   }
};

struct Less_Face : public std::binary_function<Face, Face, bool> {
   bool operator()(const Face &f1, const Face &f2) const
   {
      if ( f1.get_sorted_vertex(0) < f2.get_sorted_vertex(0) ) return true;
      if ( f1.get_sorted_vertex(0) > f2.get_sorted_vertex(0) ) return false;
      if ( f1.get_sorted_vertex(1) < f2.get_sorted_vertex(1) ) return true;
      if ( f1.get_sorted_vertex(1) > f2.get_sorted_vertex(1) ) return false;
      if ( f1.get_sorted_vertex(2) < f2.get_sorted_vertex(2) ) return true;
      if ( f1.get_sorted_vertex(2) > f2.get_sorted_vertex(2) ) return false;
      if ( f1.get_sorted_vertex(3) < f2.get_sorted_vertex(3) ) return true;
      return false;
   }
};

struct Hash_Face : public std::unary_function<Face, size_t> {
   size_t operator()(const Face &f) const
   {
      int v[4];
      f.get_ordered_vertices(v);
      return HashFNV1a<sizeof(int[4])>::eval(v);
   }
};

}  // namespace BlkC

#endif
