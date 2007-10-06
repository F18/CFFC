#ifndef _BLKC_VERTEX_INCLUDED
#define _BLKC_VERTEX_INCLUDED

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

#include <cmath>
#include <functional>

#include "BlkC_parameters.h"
#include "hash.h"

namespace BlkC
{


/*******************************************************************************
 *
 * Class Vertex
 *
 * Purpose
 * =======
 *
 *   A comparable representation of a vertex.
 *
 ******************************************************************************/

class Vertex
{
  private:
   double _x, _y, _z;

  public:
   Vertex()
      :
      _x(0.),
      _y(0.),
      _z(0.)
   { }
   Vertex(const double x, const double y, const double z)
      :
      _x(x),
      _y(y),
      _z(z)
   { }
   Vertex(const double *const p)
      :
      _x(p[0]),
      _y(p[1]),
      _z(p[2])
   { }
   // Use sythesized copy, assignment, and destructor

   double x() const { return _x; }
   double y() const { return _y; }
   double z() const { return _z; }
   void set_coord(const double x, const double y, const double z)
   {
      _x = x;
      _y = y;
      _z = z;
   }
   void set_coord(const double *const p)
   {
      _x = p[0];
      _y = p[1];
      _z = p[2];
   }
};

//--The following function objects compare the location of the vertices.
//--Equal, Less, and a Hash are defined.

struct Equal_Vertex : public std::binary_function<Vertex, Vertex, bool>
{
   bool operator()(const Vertex &v1, const Vertex &v2) const
   {
      return ( equal(v1.x(), v2.x()) &&
               equal(v1.y(), v2.y()) &&
               equal(v1.z(), v2.z()) );
   }

  private:
   bool equal(const double d1, const double d2) const
   {
      return std::fabs(d1 - d2) <
         (std::min(std::fabs(d1), std::fabs(d2)) + toler)*toler;
   }
};

struct Less_Vertex : public std::binary_function<Vertex, Vertex, bool>
{
   bool operator()(const Vertex &v1, const Vertex &v2)
   {
      if ( less_than(v1.x(), v2.x()) ) return true;
      if ( greater_than() ) return false;
      if ( less_than(v1.y(), v2.y()) ) return true;
      if ( greater_than() ) return false;
      if ( less_than(v1.z(), v2.z()) ) return true;
      return false;
   }

  private:
   double _toler;
   double _diff;
   bool less_than(const double d1, const double d2)
   {
      _diff = d1 - d2;
      _toler = (std::min(std::fabs(d1), std::fabs(d2)) + halftoler)*halftoler;
      return _toler < -_diff;
   }
   bool greater_than() const  // less_than must be called before greater_than
   {
      return _toler < _diff;
   }
};

struct Hash_Vertex : public std::unary_function<Vertex, size_t>
{
   size_t operator()(const Vertex &vertex) const
   {
      double v[3];
      v[0] = vertex.x();
      v[1] = vertex.y();
      v[2] = vertex.z();
      return HashFNV1a<sizeof(double[3])>::eval(v);
   }
};

}  // namespace BlkC

#endif
