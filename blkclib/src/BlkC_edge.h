#ifndef _BLKC_EDGE_INCLUDED
#define _BLKC_EDGE_INCLUDED

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

namespace BlkC
{


/*******************************************************************************
 *
 * Class Edge
 *
 * Purpose
 * =======
 *
 *   A comparable representation of a block edge.
 *
 ******************************************************************************/

class Edge
{
  private:
   int _v[2];
   int _si[2];                           // sorted indices to _v

  public:
   Edge(const int v0, const int v1)
   {
      _v[0] = v0; _v[1] = v1;
      sort_vertex();
   }
   Edge(const int *const v)
   {
      _v[0] = v[0]; _v[1] = v[1];
      sort_vertex();
   }
   // Use sythesized copy, assignment, and destructor

   unsigned int get_num_vertices() const { return 2; }
   int get_vertex(const int i) const { return _v[i]; }
   int get_sorted_vertex(const int i) const { return _v[_si[i]]; }
   int get_min_vertex() const { return _v[_si[0]]; }
   int get_max_vertex() const { return _v[_si[1]]; }

  private:
   void sort_vertex()
   {
      if ( _v[1] < _v[0] ) {
         _si[0] = 1;
         _si[1] = 0;
      }
      else {
         _si[0] = 0;
         _si[1] = 1;
      }
   }
};

//--Operators for comparing edge

bool operator==(const Edge &e1, const Edge &e2)
{
   return ( e1.get_min_vertex() == e2.get_min_vertex() &&
            e1.get_max_vertex() == e2.get_max_vertex() );
}

//--The following function objects compare the addresses of the vertices.
//--Equal, Less, and a Hash are defined.

struct Equal_Edge : public std::binary_function<Edge, Edge, bool> {
   bool operator()(const Edge &e1, const Edge &e2) const
   {
      return ( e1 == e2 );
   }
};

struct Less_Edge : public std::binary_function<Edge, Edge, bool> {
   bool operator()(const Edge &e1, const Edge &e2) const
   {
      if ( e1.get_min_vertex() < e2.get_min_vertex() ) return true;
      if ( e1.get_min_vertex() > e2.get_min_vertex() ) return false;
      if ( e1.get_max_vertex() < e2.get_max_vertex() ) return true;
      return false;
   }
};

struct Hash_Edge : public std::unary_function<Edge, size_t> {
   size_t operator()(const Edge &e) const
   {
      int v[2];
      v[0] = e.get_min_vertex();
      v[1] = e.get_max_vertex();
      return HashFNV1a<sizeof(int[2])>::eval(v);
   }
};

}  // namespace BlkC

#endif
