#ifndef _MESH_NEIGHBOUR_INCLUDED
#define _MESH_NEIGHBOUR_INCLUDED

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

   Adapted from GMSH source code - Copyright (C) 2006 S. Guzik, C. Geuzaine,
   J.-F. Remacle and licensed under GNU General Public License version 2.

*******************************************************************************/


/*******************************************************************************
 *
 * - The classes in this file construct a database of the blocks that share
 *   the lower-dimensional bounding objects (e.g., vertex, edge, or face).
 *   These lower-dimensional objects are referred to in general as polytopes.
 *
 ******************************************************************************/

#include <algorithm>
#include <iterator>
#include <list>
#include <map>
#include <vector>

// #define HAVE_HASH_MAP

#if defined(HAVE_HASH_MAP)
#include "hash_map.h"
#endif


/*==============================================================================
 * File scope types
 *============================================================================*/

typedef std::list<BlkC::Block*> Neighbours;
typedef Neighbours::const_iterator NeighboursConstIterator;
typedef Neighbours::iterator NeighboursIterator;

struct Range_t
{
   int num;
   NeighboursIterator begin;
   Range_t() : num(0) { }
};

//--Use a hash map for neighbour lookup if possible, otherwise a map will do

#if defined(HAVE_HASH_MAP)
struct Hash_Vertex : public std::unary_function<int, size_t>
{
   size_t operator()(const int *const v) const
   {
      return HashFNV1a<sizeof(int)>::eval(v);
   }
};
typedef HASH_MAP<int, Range_t, Hash_Vertex, std::equal_to<int> > VertNRange;
typedef HASH_MAP<BlkC::Edge, Range_t, BlkC::Hash_Edge, BlkC::Equal_Edge>
   EdgeNRange;
typedef HASH_MAP<BlkC::Block, Range_t, BlkC::Hash_Face, BlkC::Equal_Face>
   FaceNRange;
#else
typedef std::map<int, Range_t, std::less<int> > VertNRange;
typedef std::map<BlkC::Edge, Range_t, BlkC::Less_Edge> EdgeNRange;
typedef std::map<BlkC::Face, Range_t, BlkC::Less_Face> FaceNRange;
#endif
typedef VertNRange::iterator VertNRangeIterator;
typedef VertNRange::const_iterator VertNRangeConstIterator;
typedef EdgeNRange::iterator EdgeNRangeIterator;
typedef EdgeNRange::const_iterator EdgeNRangeConstIterator;
typedef FaceNRange::iterator FaceNRangeIterator;
typedef FaceNRange::const_iterator FaceNRangeConstIterator;


/*==============================================================================
 * Traits classes - that return information about a type
 *============================================================================*/

//--This is a traits/policy class for the lower-dimension polytopes that bound
//--an block (i.e., vertex, edge, or face).  It returns the corresponding
//--range type.

template <typename Polytope> struct PolytopeTr;
template <> struct PolytopeTr<int>
{
   typedef VertNRange PolytopeNRange;
   typedef VertNRangeConstIterator PNRConstIterator;
   static int getPolytope(const BlkC::Block *const block, const int nPolytope)
   {
      return block->get_vertex(nPolytope);
   }
};
template <> struct PolytopeTr<BlkC::Edge>
{
   typedef EdgeNRange PolytopeNRange;
   typedef EdgeNRangeConstIterator PNRConstIterator;
   static BlkC::Edge getPolytope(const BlkC::Block *const block,
                                 const int nPolytope)
   {
      return block->get_edge(nPolytope);
   }
};
template <> struct PolytopeTr<BlkC::Face>
{
   typedef FaceNRange PolytopeNRange;
   typedef FaceNRangeConstIterator PNRConstIterator;
   static BlkC::Face getPolytope(const BlkC::Block *const block,
                                 const int nPolytope)
   {
      return block->get_face(nPolytope);
   }
};


/*******************************************************************************
 *
 * class: PolytopeIterator
 *
 * Purpose
 * =======
 *
 *   Provides an iterator to iterate through the unique vertices, edges and,
 *   faces defined by the (hash_)maps.  An iterator to the (hash_)map has a
 *   value type of a (key, value) pair.  This iterator makes the container look
 *   like it only contains the bounding polytope.
 *
 * Constructors
 * ============
 *
 *   The class takes one template argument <Polytope> which should be either
 *   int, BlkC::Edge, or BlkC::Face.
 *
 *   PolytopeIterator()
 *                     -- default constructor
 *
 *   PolytopeIterator(const BaseIterator &baseIter_in)
 *                     -- the base iterator is the real iterator to the
 *                        (hash_)map
 *
 *   PolytopeIterator(const  PolytopeIterator &iter)
 *                     -- copy
 *
 * Destructor
 * ==========
 *
 *   ~PolytopeIterator -- synthesized
 *
 * Member Functions
 * ================
 *
 *   int num_neighbours()
 *                     -- returns number of blocks sharing the polytope
 *
 * Operators
 * =========
 *
 *   Includes typical bidirectional iterator operators {=, ==, !=, *, ->, ++(),
 *   --(), ()++, ()--} with * and -> dereferencing to the polytope.
 *
 * Notes
 * =====
 *
 *   - Only constant iterators are supported.
 *   
 ******************************************************************************/

template<typename Polytope>
class PolytopeIterator  // Actually only a const_iterator
{

  public:

   // The base iterator iterates through the (hash_)maps defining the range of
   // block neighbours for a polytope
   typedef typename PolytopeTr<Polytope>::PNRConstIterator BaseIterator;
   typedef PolytopeIterator Self;

   // Standard traits
   typedef std::bidirectional_iterator_tag iterator_category;
   typedef Polytope value_type;
   typedef BaseIterator difference_type;
   typedef const Polytope* pointer;
   typedef const Polytope& reference;

//--Constructors

   PolytopeIterator() : baseIter() { }
   PolytopeIterator(const BaseIterator &baseIter_in) : baseIter(baseIter_in) { }
   PolytopeIterator(const PolytopeIterator &iter) : baseIter(iter.baseIter) { }
   PolytopeIterator& operator=(const PolytopeIterator& iter)
   {
      if ( &iter != this ) baseIter = iter.baseIter;
      return *this;
   }

//--Comparison

   bool operator==(const Self& iter) const { return baseIter == iter.baseIter; }
   bool operator!=(const Self& iter) const { return baseIter != iter.baseIter; }

//--Dereference and arrow operator

   reference operator*() const { return baseIter->first; }
   pointer operator->() const { return &baseIter->first; }

//--Increment

   // Prefix
   Self& operator++()
   {
      ++baseIter;
      return *this;
   }

   // Postfix
   Self operator++(int)
   {
      Self tmp = *this;
      ++*this;
      return tmp;
   }

//--Decrement

   // Prefix
   Self& operator--()
   {
      --baseIter;
      return *this;
   }

   // Postfix
   Self operator--(int)
   {
      Self tmp = *this;
      --*this;
      return tmp;
   }

//--Special

   int num_neighbours() const { return baseIter->second.num; }

  private:
   NeighboursConstIterator begin_neighbours() const
   {
      return baseIter->second.begin;
   }

//--Data members

   BaseIterator baseIter;

//--Friends

   friend class MNeighbour;

};


/*******************************************************************************
 *
 * Template meta-programming classes
 *
 ******************************************************************************/

/*==============================================================================
 *
 * struct FindNeighbours
 *
 * Purpose
 * =======
 *
 *   This is the main working routine for finding the blocks sharing a lower-
 *   dimension polytope (i.e, the structures bounding the block: vertex, edge,
 *   or face).  It is templated for any type of polytope and uses template meta-
 *   programming to unroll the loop over the number of that type of polytope in
 *   the block.
 *
 *============================================================================*/

//--Entry point

template <typename Polytope,            // Polytope is a lower dimensional
                                        // bounding object (int, Edge, or
                                        // Face) of the block.  We want to
                                        // find the blocks sharing a given
                                        // Polytope.
          int NP>                       // NP is an index through the number
                                        // of Polytope in the block (e.g.,
                                        // number of vertices)
struct FindNeighbours
{
   typedef typename PolytopeTr<Polytope>::PolytopeNRange PolytopeNRange;
   static void eval(BlkC::Block *block, PolytopeNRange &polytopeNRange,
                    Neighbours &neighbours)
   {
      Polytope polytope = PolytopeTr<Polytope>::getPolytope(block, NP-1);
      Range_t &range = polytopeNRange[polytope];
      if ( range.num == 0 ) range.begin = neighbours.end();
      range.begin = neighbours.insert(range.begin, block);
      ++range.num;
      FindNeighbours<Polytope, NP-1>::eval(block, polytopeNRange, neighbours);
   }
};

//--Terminate loop when NP == 1

template <typename Polytope>
struct FindNeighbours<Polytope, 1>
{
   typedef typename PolytopeTr<Polytope>::PolytopeNRange PolytopeNRange;
   static void eval(BlkC::Block *block, PolytopeNRange &polytopeNRange,
                    Neighbours &neighbours)
   {
      Polytope polytope = PolytopeTr<Polytope>::getPolytope(block, 0);
      Range_t &range = polytopeNRange[polytope];
      if ( range.num == 0 ) range.begin = neighbours.end();
      range.begin = neighbours.insert(range.begin, block);
      ++range.num;
   }
};


/*******************************************************************************
 *
 * class: MNeighbour
 *
 * Purpose
 * =======
 *
 *   Generates a database for quick lookup of blocks sharing a vertex (int),
 *   Edge, or Face (referred to in general as polytopes.  I.e., a geometrical
 *   construct bounding the block with a lower dimension than the block).  The
 *   neighbour blocks are stored in a list.  Maps or hash maps provide an index
 *   into the list for a given polytope.  In addition to the block neighbours,
 *   this class provides iterators to the unique vertices, edges, and faces in
 *   the database.
 *
 * Iterators
 * =========
 *
 *   Vertex_const_iterator
 *                     -- iterator to sequence of unique vertices
 *   Edge_const_iterator
 *                     -- iterator to sequence of unique edges
 *   Face_const_iterator
 *                     -- iterator to sequence of unique faces
 *
 * Constructors
 * ============
 *
 *   MNeighbour()      -- default constructor does nothing
 *
 * Destructor
 * ==========
 *
 *   ~MNeighbour()     -- synthesized
 *
 * Member Functions
 * ================
 *
 * Iterators
 * ---------
 *
 *   Vertex_const_iterator vertex_begin()
 *                     -- first unique vertex
 *
 *   Vertex_const_iterator vertex_end()
 *                     -- one-past-last unique vertex
 *
 *   Edeg_const_iterator edge_begin()
 *                     -- first unique edge
 *
 *   Edge_const_iterator edge_end()
 *                     -- one-past-last unique edge
 *
 *   Face_const_iterator face_begin()
 *                     -- first unique face
 *
 *   Face_const_iterator face_end()
 *                     -- one-past-last unique face
 *
 * Routines for adding blocks
 * --------------------------
 *
 *   void add_block(BlkC::Block *const block)
 *                     -- adds a block to the neighbours database.
 *     block              (I) representation of the block
 *
 *   void clear()      -- clears the neighbours database.  Use this before
 *                        constructing a new database of block neighbours.
 *                        Otherwise, new blocks are added to the existing
 *                        database.
 *
 * Routines for querying the neighbours database
 * ---------------------------------------------
 *
 *   int max_vertex_neighbours()
 *                     -- returns max number of blocks sharing a vertex.  This
 *                        is also the max number of blocks sharing any
 *                        bounding polytope.
 *
 *   int max_edge_neighbours()
 *                     -- returns max number of blocks sharing an edge.
 *
 *   int max_face_neighbours()
 *                     -- returns max number of blocks sharing a face.
 *
 *   int vertex_iterator(const int vertex, Vertex_const_iterator &itVert)
 *                     -- find the polytope iterator for a vertex
 *     vertex             (I) vertex to find
 *     itVert             (O) iterator to the vertex
 *     return             0 - success
 *                        1 - not found
 *
 *   int edge_iterator(const BlkC::Edge &edge, Edge_const_iterator &itEdge)
 *                     -- find the polytope iterator for an edge
 *     edge               (I) edge to find
 *     itVert             (O) iterator to the edge
 *     return             0 - success
 *                        1 - not found
 *
 *   int face_iterator(const BlkC::Face &face, Face_const_iterator &itFace)
 *                     -- find the polytope iterator for a face
 *     face               (I) face to find
 *     itVert             (O) iterator to the face
 *     return             0 - success
 *                        1 - not found
 *
 *   int vertex_neighbours(Vertex_const_iterator &itVert, BlkC::Block **pBlk)
 *                     -- return all blocks sharing a vertex
 *     itVert             (I) iterator to a unique vertex.  Avoids a find()
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks
 *
 *   int vertex_neighbours(const int vertex, BlkC::Block **pBlk)
 *                     -- returns all blocks sharing a vertex
 *     vertex             (I) vertex.  Must be found in (hash_)map
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks.  0 if vertex not found
 *
 *   int edge_neighbours(Edge_const_iterator &itEdge, BlkC::Block **pBlk)
 *                     -- return all blocks sharing a edge
 *     itEdge             (I) iterator to a unique edge.  Avoids a find()
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks
 *
 *   int edge_neighbours(const BlkC::Edge &edge, BlkC::Block **pBlk)
 *                     -- returns all blocks sharing a edge
 *     edge               (I) edge.  It must be found in (hash_)map
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks.  0 if edge not found
 *
 *   int face_neighbours(Face_const_iterator &itFace, BlkC::Block **pBlk)
 *                     -- return all blocks sharing a face
 *     itFace             (I) iterator to a unique face.  Avoids a find()
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks
 *
 *   int face_neighbours(const BlkC::Face &face, BlkC::Block **pBlk)
 *                     -- returns all blocks sharing a face
 *     face               (I) face.  It must be found in (hash_)map
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of neighbour blocks.  0 if face not found
 *
 *   int vertex_num_neighbours(const int vertex)
 *                     -- returns the number of blocks sharing a vertex
 *
 *   int edge_num_neighbours(const BlkC::Edge& edge)
 *                     -- returns the number of blocks sharing an edge
 *
 *   int face_num_neighbours(const BlkC::Face& face)
 *                     -- returns the number of blocks sharing a face
 *
 *   int all_block_neighbours(BlkC::Block *const block,
 *                            BlkC::Block **pBlk
 *                            const bool exclusive = true)
 *                     -- find all the blocks connected to the given block
 *     block              (I) block to find neighbours to
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     exclusive          (I) T - exclude the given block in vBlk
 *                            F - include the given block in vBlk
 *     return             number of blocks added to pBlk
 *
 *   int all_block_edge_neighbours(BlkC::Block *const block,
 *                                 BlkC::Block **pBlk
 *                                 const bool exclusive = true)
 *                     -- find all the blocks connected to the edges of a
 *                        given block
 *     block              (I) block to find neighbours to
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     exclusive          (I) T - exclude the given block in vBlk
 *                            F - include the given block in vBlk
 *     return             number of blocks added to pBlk
 *
 *   int all_block_face_neighbours(BlkC::Block *const block,
 *                                 BlkC::Block **pBlk
 *                                 const bool exclusive = true)
 *                     -- find all the blocks connected to the faces of a
 *                        given block
 *     block              (I) block to find neighbours to
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     exclusive          (I) T - exclude the given block in vBlk
 *                            F - include the given block in vBlk
 *     return             number of blocks added to pBlk
 *
 *   int block_vertex_neighbours(const BlkC::Block *const block,
 *                               const Vertex_const_iterator &itVert,
 *                               BlkC::Block **pBlk)
 *                     -- find all the blocks connected to a vertex of a given
 *                        block.  This is the same as routine
 *                        vertex_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to exclude
 *     itVert             (I) iterator to a unique vertex.  Avoids a find()
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of blocks added to pBlk
 *
 *   int block_vertex_neighbours(const BlkC::Block *const block,
 *                               const int vertex, BlkC::Block **pBlk)
 *                     -- find all the blocks connected to a vertex of a given
 *                        block.  This is the same as routine
 *                        vertex_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to exclude
 *     vertex             (I) vertex to find block neighbours to
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of blocks added to pBlk
 *
 *   int block_edge_neighbours(const BlkC::Block *const block,
 *                             const Edge_const_iterator &itEdge,
 *                             BlkC::Block **pBlk)
 *                     -- find all the blocks connected to an edge of a given
 *                        block.  This is the same as routine
 *                        edge_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to exclude
 *     itEdge             (I) iterator to a unique edge.  Avoids a find()
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of blocks added to pBlk
 *
 *   int block_edge_neighbours(const BlkC::Block *const block,
 *                             const BlkC::Edge &edge, BlkC::Block **pBlk)
 *                     -- find all the blocks connected to an edge of a given
 *                        block.  This is the same as routine
 *                        edge_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to exclude
 *     edge               (I) edge to find block neighbours to
 *     pBlk               (O) array with blocks loaded starting at pBlk[0]
 *     return             number of blocks added to pBlk
 *
 *   BlkC::Block *block_face_neighbour(const BlkC::Block *const block,
 *                                     const Face_const_iterator &itFace)
 *                     -- find the block connected to the face of a given
 *                        block.  This is the same as routine
 *                        face_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to find neighbour to
 *     itFace             (I) iterator to a unique face.  Avoids a find()
 *     return             the neighbour block or 0 if none found
 *
 *   BlkC::Block *block_face_neighbour(const BlkC::Block *const block,
 *                                     const BlkC::Face &face)
 *                     -- find the block connected to the face of a given
 *                        block.  This is the same as routine
 *                        face_neighbours except that the given block is
 *                        excluded from the result.
 *     block              (I) block to find neighbour to
 *     face               (I) face to find block neighbour across
 *     return             the neighbour block or 0 if none found
 *
 ******************************************************************************/

class MNeighbour
{


/*==============================================================================
 * Class scope types
 *============================================================================*/

  private:
   typedef std::set<BlkC::Block*, std::less<BlkC::Block*> > BlkSet;
   typedef BlkSet::const_iterator BlkSetConstIterator;


/*==============================================================================
 * Iterators
 *============================================================================*/

  public:

   // Iterators over polytopes
   typedef PolytopeIterator<int> Vertex_const_iterator;
   typedef PolytopeIterator<BlkC::Edge> Edge_const_iterator;
   typedef PolytopeIterator<BlkC::Face> Face_const_iterator;

   // Begin and end points for these iterators
   Vertex_const_iterator vertex_begin() const
   {
      return Vertex_const_iterator(vertNRange.begin()); 
   }
   Vertex_const_iterator vertex_end() const
   {
      return Vertex_const_iterator(vertNRange.end()); 
   }
   Edge_const_iterator edge_begin() const
   {
      return Edge_const_iterator(edgeNRange.begin()); 
   }
   Edge_const_iterator edge_end() const
   {
      return Edge_const_iterator(edgeNRange.end()); 
   }
   Face_const_iterator face_begin() const
   {
      return Face_const_iterator(faceNRange.begin()); 
   }
   Face_const_iterator face_end() const
   {
      return Face_const_iterator(faceNRange.end()); 
   }


/*==============================================================================
 * Public member functions and data
 *============================================================================*/

  public:

//--Default constructor.

   // The constructor cannot really do anything because the initialization
   // functions (find_neighbours_*) cannot perform template argument deduction.
   MNeighbour()
      :
      maxVertNeighbours(0),
      maxEdgeNeighbours(0),
      maxFaceNeighbours(0)
   { }

//--Add a block to the database

   void add_block(BlkC::Block *const block)
   {
      // Find neighbours for vertices
      FindNeighbours<int, 8>::eval(block, vertNRange, neighbours);
      // Find neighbours for edges
      FindNeighbours<BlkC::Edge, 12>::eval(block, edgeNRange, neighbours);
      // Find neighbours for faces
      FindNeighbours<BlkC::Face, 6>::eval(block, faceNRange, neighbours);
      // Flag that max neighbours must be computed
      maxVertNeighbours = -1;
      maxEdgeNeighbours = -1;
      maxFaceNeighbours = -1;
   }

/*--------------------------------------------------------------------*
 * Reset the database
 *--------------------------------------------------------------------*/

   void clear() 
   {
      vertNRange.clear();
      edgeNRange.clear();
      faceNRange.clear();
      neighbours.clear();
      maxVertNeighbours = 0;
      maxEdgeNeighbours = 0;
      maxFaceNeighbours = 0;
   }

/*--------------------------------------------------------------------*
 * Query the database
 *--------------------------------------------------------------------*/

//--Get the max number of blocks sharing a vertex, edge, or face

   int max_vertex_neighbours()
   {
      if ( maxVertNeighbours < 0 ) compute_max_vert_neighbours();
      return maxVertNeighbours;
   }
   int max_edge_neighbours()
   {
      if ( maxEdgeNeighbours < 0 ) compute_max_edge_neighbours();
      return maxEdgeNeighbours;
   }
   int max_face_neighbours()
   {
      if ( maxFaceNeighbours < 0 ) compute_max_face_neighbours();
      return maxFaceNeighbours;
   }

//--Return iterators to polytopes

   // Vertex
   int vertex_iterator(const int vertex, Vertex_const_iterator &itVert) const
   {
      VertNRangeConstIterator itVNR = vertNRange.find(vertex);
      if ( itVNR == vertNRange.end() ) return 1;  // Vertex not found
      itVert.baseIter = itVNR;
      return 0;
   }

   // Edge
   int edge_iterator(const BlkC::Edge &edge, Edge_const_iterator &itEdge) const
   {
      EdgeNRangeConstIterator itENR = edgeNRange.find(edge);
      if ( itENR == edgeNRange.end() ) return 1;  // Edge not found
      itEdge.baseIter = itENR;
      return 0;
   }

   // Face
   int face_iterator(const BlkC::Face &face, Face_const_iterator &itFace) const
   {
      FaceNRangeConstIterator itFNR = faceNRange.find(face);
      if ( itFNR == faceNRange.end() ) return 1;  // Face not found
      itFace.baseIter = itFNR;
      return 0;
   }

//--Return blocks that share a vertex

   // From a vertex iterator
   int vertex_neighbours(Vertex_const_iterator &itVert, BlkC::Block **pBlk)
      const
   {
      NeighboursConstIterator itBlk = itVert.begin_neighbours();
      for ( int n = itVert.num_neighbours(); n--; ) *pBlk++ = *itBlk++;
      return itVert.num_neighbours();
   }

   // From a vertex
   int vertex_neighbours(const int vertex, BlkC::Block **pBlk) const
   {
      VertNRangeConstIterator itVert = vertNRange.find(vertex);
      if ( itVert == vertNRange.end() ) return 0;  // Vertex not found
      NeighboursConstIterator itBlk = itVert->second.begin;
      for ( int n = itVert->second.num; n--; ) *pBlk++ = *itBlk++;
      return itVert->second.num;
   }

//--Return blocks that share an edge

   // From an edge iterator
   int edge_neighbours(Edge_const_iterator &itEdge, BlkC::Block **pBlk) const
   {
      NeighboursConstIterator itBlk = itEdge.begin_neighbours();
      for ( int n = itEdge.num_neighbours(); n--; ) *pBlk++ = *itBlk++;
      return itEdge.num_neighbours();
   }

   // From an edge
   int edge_neighbours(const BlkC::Edge &edge, BlkC::Block **pBlk) const
   {
      EdgeNRangeConstIterator itEdge = edgeNRange.find(edge);
      if ( itEdge == edgeNRange.end() ) return 0;  // Edge not found
      NeighboursConstIterator itBlk = itEdge->second.begin;
      for ( int n = itEdge->second.num; n--; ) *pBlk++ = *itBlk++;
      return itEdge->second.num;
   }

//--Return blocks that share a face

   // From a face iterator
   int face_neighbours(Face_const_iterator &itFace, BlkC::Block **pBlk) const
   {
      NeighboursConstIterator itBlk = itFace.begin_neighbours();
      for ( int n = itFace.num_neighbours(); n--; )
         *pBlk++ = *itBlk++;
      return itFace.num_neighbours();
   }

   // From a face
   int face_neighbours(const BlkC::Face &face, BlkC::Block **pBlk) const
   {
      FaceNRangeConstIterator itFace = faceNRange.find(face);
      if ( itFace == faceNRange.end() ) return 0;  // Face not found
      NeighboursConstIterator itBlk = itFace->second.begin;
      for ( int n = itFace->second.num; n--; ) *pBlk++ = *itBlk++;
      return itFace->second.num;
   }

//--Return the number of blocks that share a vertex

   int vertex_num_neighbours(const int vertex) const
   {
      VertNRangeConstIterator itVert = vertNRange.find(vertex);
      if ( itVert == vertNRange.end() ) return 0;  // Vertex not found
      return itVert->second.num;
   }

//--Return the number of blocks that share an edge

   int edge_num_neighbours(const BlkC::Edge& edge) const
   {
      EdgeNRangeConstIterator itEdge = edgeNRange.find(edge);
      if ( itEdge == edgeNRange.end() ) return 0;  // Edge not found
      return itEdge->second.num;
   }

//--Return the number of blocks that share a face

   int face_num_neighbours(const BlkC::Face& face) const
   {
      FaceNRangeConstIterator itFace = faceNRange.find(face);
      if ( itFace == faceNRange.end() ) return 0;  // Face not found
      return itFace->second.num;
   }

//--Return all the blocks connected to a given block

   int all_block_neighbours(BlkC::Block *const block,
                            BlkC::Block **pBlk,
                            const bool exclusive = true)
   {
      block_neighbours1<int>(block, vertNRange, 8);
      if ( exclusive ) blkSet.erase(block);
      const BlkSetConstIterator itBlkEnd = blkSet.end();
      for ( BlkSetConstIterator itBlk = blkSet.begin(); itBlk != itBlkEnd;
            ++itBlk ) *pBlk++ = *itBlk;
      const int nBlk = blkSet.size();
      blkSet.clear();
      return nBlk;
   }

//--Return all the blocks connected to the edges of a given block

   int all_block_edge_neighbours(BlkC::Block *const block,
                                 BlkC::Block **pBlk,
                                 const bool exclusive = true)
   {
      block_neighbours1<BlkC::Edge>(block, edgeNRange, 12);
      if ( exclusive ) blkSet.erase(block);
      const BlkSetConstIterator itBlkEnd = blkSet.end();
      for ( BlkSetConstIterator itBlk = blkSet.begin(); itBlk != itBlkEnd;
            ++itBlk ) *pBlk++ = *itBlk;
      const int nBlk = blkSet.size();
      blkSet.clear();
      return nBlk;
   }

//--Return all the blocks connected to the faces of a given block

   int all_block_face_neighbours(BlkC::Block *const block,
                                 BlkC::Block **pBlk,
                                 const bool exclusive = true)
   {
      block_neighbours1<BlkC::Face>(block, faceNRange, 6);
      if ( exclusive ) blkSet.erase(block);
      const BlkSetConstIterator itBlkEnd = blkSet.end();
      for ( BlkSetConstIterator itBlk = blkSet.begin(); itBlk != itBlkEnd;
            ++itBlk ) *pBlk++ = *itBlk;
      const int nBlk = blkSet.size();
      blkSet.clear();
      return nBlk;
   }

//--Return blocks connected to a specific vertex of a given block (Note:
//--these routines are the same as vertex_neighbours() except that the given
//--block is excluded)

   // From a vertex interator
   int block_vertex_neighbours(const BlkC::Block *block,
                               const Vertex_const_iterator &itVert,
                               BlkC::Block **pBlk) const
   {
      NeighboursConstIterator itBlk = itVert.begin_neighbours();
      int n = itVert.num_neighbours();
      while ( *itBlk != block ) {
         *pBlk++ = *itBlk++;
         --n;
      }
      ++itBlk;
      for ( --n; n--; ) *pBlk++ = *itBlk++;
      return itVert.num_neighbours() - 1;
   }

   // From a vertex
   int block_vertex_neighbours(const BlkC::Block *block,
                               const int vertex,
                               BlkC::Block **pBlk) const
   {
      VertNRangeConstIterator itVert = vertNRange.find(vertex);
      if ( itVert == vertNRange.end() ) return 0;  // Vertex not found
      return
         block_vertex_neighbours(block, Vertex_const_iterator(itVert), pBlk);
   }

//--Return blocks connected to a specific edge of a given block (Note: these
//--routines are the same as edge_neighbours except that the given block is
//--excluded)

   // From an edge iterator
   int block_edge_neighbours(const BlkC::Block *const block,
                             const Edge_const_iterator &itEdge,
                             BlkC::Block **pBlk) const
   {
      NeighboursConstIterator itBlk = itEdge.begin_neighbours();
      int n = itEdge.num_neighbours();
      while( *itBlk != block ) {
         *pBlk++ = *itBlk++;
         --n;
      }
      ++itBlk;
      for ( --n; n--; ) *pBlk++ = *itBlk++;
      return itEdge.num_neighbours() - 1;
   }

   // From an edge
   int block_edge_neighbours(const BlkC::Block *const block,
                             const BlkC::Edge &edge,
                             BlkC::Block **pBlk) const
   {
      EdgeNRangeConstIterator itEdge = edgeNRange.find(edge);
      if ( itEdge == edgeNRange.end() ) return 0;  // Edge not found
      return block_edge_neighbours(block,  Edge_const_iterator(itEdge), pBlk);
   }

//--Return the block connected to a specific face of a given block (Note: these
//--routines are the same as face_neighbours except that the given block is
//--excluded)

   // From a face iterator
   BlkC::Block *block_face_neighbour(const BlkC::Block *const block,
                                     const Face_const_iterator &itFace) const
   {
      NeighboursConstIterator itBlk = itFace.begin_neighbours();
      if ( *itBlk != block ) return *itBlk;
      if ( itFace.num_neighbours() > 1 ) {
         ++itBlk;
         if ( *itBlk != block ) return *itBlk;
      }
      return 0;
   }

   // From a face
   BlkC::Block *block_face_neighbour(const BlkC::Block *const block,
                                     const BlkC::Face &face) const
   {
      FaceNRangeConstIterator itFace = faceNRange.find(face);
      if ( itFace == faceNRange.end() ) return 0;  // Face not found
      return block_face_neighbour(block, Face_const_iterator(itFace));
   }


/*==============================================================================
 * Private member functions and data
 *============================================================================*/

private:

//--Data members

   VertNRange vertNRange;               // Map of Vert. to indexes in neighbours
   EdgeNRange edgeNRange;               // Map of edge to indexes in neighbours
   FaceNRange faceNRange;               // Map of face to indexes in neighbours
   Neighbours neighbours;               // List of block neighbours
   BlkSet blkSet;                     // Buffer for finding unique blocks
   int maxVertNeighbours;
   int maxEdgeNeighbours;
   int maxFaceNeighbours;

//--Compute max neighbours

   void compute_max_vert_neighbours()
   {
      maxVertNeighbours = 0;
      for ( VertNRangeConstIterator itVert = vertNRange.begin();
            itVert != vertNRange.end(); ++itVert )
         maxVertNeighbours = std::max(maxVertNeighbours, itVert->second.num);
   }
   void compute_max_edge_neighbours()
   {
      maxEdgeNeighbours = 0;
      for ( EdgeNRangeConstIterator itEdge = edgeNRange.begin();
            itEdge != edgeNRange.end(); ++itEdge )
         maxEdgeNeighbours = std::max(maxEdgeNeighbours, itEdge->second.num);
   }
   void compute_max_face_neighbours()
   {
      maxFaceNeighbours = 0;
      for ( FaceNRangeConstIterator itFace = faceNRange.begin();
            itFace != faceNRange.end(); ++itFace )
         maxFaceNeighbours = std::max(maxFaceNeighbours, itFace->second.num);
   }

//--Work routine for finding all the block neighbours

   template <typename Polytope>
   void block_neighbours1(const BlkC::Block *const block,
                          typename PolytopeTr<Polytope>::PolytopeNRange
                          &polytopeNRange,
                          const int nPPE)
   {
      for ( int iPPE = 0; iPPE != nPPE; ++iPPE ) {
         Polytope polytope = PolytopeTr<Polytope>::getPolytope(block, iPPE);
         typename PolytopeTr<Polytope>::PNRConstIterator itPoly =
            polytopeNRange.find(polytope);
         if ( itPoly != polytopeNRange.end() ) {
            NeighboursConstIterator itBlk = itPoly->second.begin;
            for ( int n = itPoly->second.num; n--; )
               blkSet.insert(*itBlk++);
         }
      }
   }

};

#endif
