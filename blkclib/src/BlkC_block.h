#ifndef _BLKC_BLOCK_INCLUDED
#define _BLKC_BLOCK_INCLUDED

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

#include <cstring>

#include "BlkC_edge.h"
#include "BlkC_face.h"

namespace BlkC
{

struct NeighbourData
{
   int label;                           // Label for the neighbour block
   int boundary[3];                     // Boundary of the neighbour block
   int ctm[6];                          // Compact transformation matrix
                                        // (including offset) for transformation
                                        // to the neighbour block
   NeighbourData *next;
   NeighbourData() : next(0) { }
};


/*******************************************************************************
 *
 * Class Block
 *
 * Purpose
 * =======
 *
 *   A representation of a block given by integer representation of corner
 *   vertices.
 *
 * Notes
 * =====
 *
 *   - Conventions for node numbering in the block
 *
 *         2----------3
 *        /|         /|
 *       / |        / |          ^ j
 *      /  |       /  |          |
 *     6---+------7   |          |
 *     |   |      |   |          |      i
 *     |   0------+---1          +------>
 *     |  /       |  /          /
 *     | /        | /          /
 *     |/         |/          / k
 *     4----------5          `-
 *
 *   - Static const data types are defined in "block_connectivity.cc"
 *
 ******************************************************************************/

class Block
{

  public:
   int label;                           // Label or index of this block
   int dim[3];                          // Dimensions of this block.  What is
                                        // stored is dim-1, i.e., the step to
                                        // go from a cell on the min face to a
                                        // cell on the max face

  private:
   NeighbourData *nbrData[3][3][3];     // Info about neighbour blocks that
                                        // connect to this block across a given
                                        // boundary
   int extent[3][3][3];                 // Info about the domain extents
   int _v[8];                           // Corner vertices of this block
   double _vec[8][3][3];                // Vectors describing each edge leaving
                                        // a corner vertex

  public:
   static const int mEdge[12][2];       // Description of edges from vertices
   static const int mFace[6][4];        // Description of faces from vertices


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  public:

   // Default constructor
   Block()
      :
      label(-1)
   {
      dim[0] = 0;
      dim[1] = 0;
      dim[2] = 0;
      NeighbourData **pnd = nbrData[0][0];
      int *pe = extent[0][0];
      for ( int n = 27; n--; ) {
         *pnd++ = 0;
         *pe++ = -1;
      }
      for ( int i = 0; i != 8; ++i ) _v[i] = 0;
   }

   // Constructor
   Block(const int lbl, const int iBlkDim, const int jBlkDim, const int kBlkDim,
         const int *const v, const double *const vec)
      :
      label(lbl)
   {
      dim[0] = iBlkDim;
      dim[1] = jBlkDim;
      dim[2] = kBlkDim;
      NeighbourData **pnd = nbrData[0][0];
      int *pe = extent[0][0];
      for ( int n = 27; n--; ) {
         *pnd++ = 0;
         *pe++ = -1;
      }
      std::memcpy(_v, v, 8*sizeof(int));
      std::memcpy(_vec, vec, 72*sizeof(double));
   }
   // Use sythesized copy, assignment, and destructor but beware that the
   // pointers in NeighbourData are simply copied


/*==============================================================================
 * Member functions
 *============================================================================*/

/*--------------------------------------------------------------------*
 * Routines for setting values
 *--------------------------------------------------------------------*/

   // Set the value for the vertices
   void set_vertex(const int mv, const int val)
   {
      _v[mv] = val;
   }
   void set_vertex(const ILoc_t iloc, const JLoc_t jloc, const KLoc_t kloc,
                   const int val)
   {
      _v[ijktomv(iloc, jloc, kloc)] = val;
   }

   // Set the values for the edge vectors
   void set_vector(const double *const vec)
   {
      std::memcpy(_vec, vec, 72*sizeof(double));
   }

   // Set the neighbour data
   void set_neighbour_data(const int iloc, const int jloc, const int kloc,
                           NeighbourData *const pNbrData)
   {
      NeighbourData *prev = 0;
      NeighbourData **p = &nbrData[iloc+1][jloc+1][kloc+1];
      while ( *p ) {
         prev = *p;
         p = &(prev->next);
      }
      *p = pNbrData;
      if ( prev ) prev->next = *p;
   }

   // Set an extent
   void set_extent(const int iloc, const int jloc, const int kloc,
                   const bool atExtent)
   {
      extent[iloc+1][jloc+1][kloc+1] = atExtent;
   }

/*--------------------------------------------------------------------*
 * Information about vertices
 *--------------------------------------------------------------------*/

   // Get the vertex at cloc
   int get_vertex(const int mv) const
   {
      return _v[mv];
   }

   // Get the vertex at iloc, jloc, kloc
   int get_vertex(const ILoc_t iloc, const JLoc_t jloc, const KLoc_t kloc) const
   {
      return _v[ijktomv(iloc, jloc, kloc)];
   }

   // Given a vertex, find the vertex index mv
   int find_vertex_index(const int vert) const
   {
      for ( int mv = 0; mv != 8; ++mv )
         if ( vert == get_vertex(mv) ) return mv;
      return -1;  // Not found
   }

/*--------------------------------------------------------------------*
 * Information about edges
 *--------------------------------------------------------------------*/

   // Get the edge me
   Edge get_edge(const int me) const
   {
      return Edge(_v[mEdge[me][0]], _v[mEdge[me][1]]);
   }

   // Get the edge at iloc, jloc, kloc
   Edge get_edge(const ILoc_t iloc, const JLoc_t jloc, const KLoc_t kloc) const
   {
      const int me = ijktome(iloc, jloc, kloc);
      return Edge(_v[mEdge[me][0]], _v[mEdge[me][1]]);
   }

   // Given an edge, find the edge index me
   int find_edge_index(const Edge& edge) const
   {
      for ( int me = 0; me != 12; ++me )
         if ( edge == get_edge(me) ) return me;
      return -1;  // Not found
   }

/*--------------------------------------------------------------------*
 * Information about faces
 *--------------------------------------------------------------------*/

   // Get the face mf
   Face get_face(const int mf) const
   {
      return Face(_v[mFace[mf][0]], _v[mFace[mf][1]],
                  _v[mFace[mf][2]], _v[mFace[mf][3]]);
   }

   // Get the face at iloc, jloc, kloc
   Face get_face(const int iloc, const int jloc, const int kloc) const
   {
      const int mf = ijktomf(iloc, jloc, kloc);
      return Face(_v[mFace[mf][0]], _v[mFace[mf][1]],
                  _v[mFace[mf][2]], _v[mFace[mf][3]]);
   }

   // Given a face, find the face index mf
   int find_face_index(const Face& face) const
   {
      for ( int mf = 0; mf != 6; ++mf )
         if ( face == get_face(mf) ) return mf;
      return -1;  // Not found
   }

/*--------------------------------------------------------------------*
 * Information about the neighbours
 *--------------------------------------------------------------------*/

   // Check if any neighbour exists at a boundary
   bool have_neighbour(const int iloc, const int jloc, const int kloc) const
   {
      return nbrData[iloc+1][jloc+1][kloc+1];
   }

   // Check if a specific neighbour exists at a boundary
   bool have_neighbour(const int iloc, const int jloc, const int kloc,
                       const int labelNbr) const
   {
      const NeighbourData *p = nbrData[iloc+1][jloc+1][kloc+1];
      while ( p ) {
         if ( p->label == labelNbr ) return true;
         p = p->next;
      }
      return false;
   }

   // Number of neighbours
   int num_neighbour(const int iloc, const int jloc, const int kloc) const
   {
      const NeighbourData *p = nbrData[iloc+1][jloc+1][kloc+1];
      int n = 0;
      while ( p ) {
         p = p->next;
         ++n;
      }
      return n;
   }

   // Neighbour data (returns 1 = fail | 0 = success)
   int neighbour_data(const int iloc, const int jloc, const int kloc,
                      int nbrloc, int &labelNbr, int &ilocNbr,
                      int &jlocNbr, int &klocNbr, int *const ctm) const
   {
      const NeighbourData *p = nbrData[iloc+1][jloc+1][kloc+1];
      while ( nbrloc-- ) {
         if ( !p ) return 1;
         p = p->next;
      }
      if ( !p ) return 1;
      labelNbr = p->label;
      ilocNbr = p->boundary[0];
      jlocNbr = p->boundary[1];
      klocNbr = p->boundary[2];
      ctm[0] = p->ctm[0];
      ctm[1] = p->ctm[1];
      ctm[2] = p->ctm[2];
      ctm[3] = p->ctm[3];
      ctm[4] = p->ctm[4];
      ctm[5] = p->ctm[5];
      return 0;
   }

/*--------------------------------------------------------------------*
 * Information about the extent
 *--------------------------------------------------------------------*/

   bool extent_known(const int iloc, const int jloc, const int kloc) const
   {
      if ( extent[iloc+1][jloc+1][kloc+1] == -1 ) return false;
      return true;
   }

   bool at_extent(const int iloc, const int jloc, const int kloc) const
   {
      return extent[iloc+1][jloc+1][kloc+1];
   }

/*--------------------------------------------------------------------*
 * Information about the vector
 *--------------------------------------------------------------------*/

   void get_vector(const int mv, const int dir, double *const vec) const
   {
      vec[0] = _vec[mv][dir][0];
      vec[1] = _vec[mv][dir][1];
      vec[2] = _vec[mv][dir][2];
   }

/*--------------------------------------------------------------------*
 * Conversions between ijk indices and a one dimensional array index
 *--------------------------------------------------------------------*/

   // Convert ijk indices to vertex index "mv" in array
   static int ijktomv(const int iloc, const int jloc, const int kloc)
   {
      // Note that + 1 >> 1 converts -1 to 0
      return
         kloc + 1 << 1 |
         jloc + 1 |
         iloc + 1 >> 1;
   }

   // Convert vertex index "mv" in array to ijk indices
   static void mvtoijk(const int cloc, int &iloc, int &jloc, int &kloc)
   {
      // Note that (<< 1) - 1 converts 0 to -1
      iloc = ((cloc & 1) << 1) - 1;
      jloc = ((cloc >> 1 & 1) << 1) - 1;
      kloc = (cloc >> 2 << 1) - 1;
   }

   // Convert ijk indices to edge index "me" in array
   static int ijktome(const int iloc, const int jloc, const int kloc)
   {
      // me gives 0, 4, or 8
      const int me = !std::abs(kloc) << 3 | !std::abs(jloc) << 2;
      return me |
         // kloc, if != 0, always shifts by 1
          kloc + 1 >> 1 << 1 |
         // jloc, if != 0, can shift by 0 or 1
         jloc + 1 >> 1 << (( iloc ) ? 1 : 0) |
         // iloc, if != 0, always shifts by 0
         iloc + 1 >> 1;
   }

   // Convert edge index "me" in array to ijk indices
   static void metoijk(const int eloc, int &iloc, int &jloc, int &kloc)
   {
      const int mm0 = ((eloc & 1) << 1) - 1;
      const int mm1 = ((eloc >> 1 & 1) << 1) - 1;
      switch ( eloc >> 2 ) {
      case 0:
         iloc = 0;
         jloc = mm0;
         kloc = mm1;
         break;
      case 1:
         iloc = mm0;
         jloc = 0;
         kloc = mm1;
         break;
      case 2:
         iloc = mm0;
         jloc = mm1;
         kloc = 0;
         break;
      }
   }

   // Convert ijk indices to face index "mf" in array
   static int ijktomf(const int iloc, const int jloc, const int kloc)
   {
      // bit defines the last bit
      const int bit = iloc + jloc + kloc;
      return
         std::abs(kloc) << 2 |
         std::abs(jloc) << 1 |
         bit + 1 >> 1;
   }

   // Convert face index "mf" in array to ijk indices
   static void mftoijk(const int floc, int &iloc, int&jloc, int &kloc)
   {
      const int mm0 = ((floc & 1) << 1) - 1;
      switch ( floc >> 1 ) {
      case 0:
         iloc = mm0;
         jloc = 0;
         kloc = 0;
         break;
      case 1:
         iloc = 0;
         jloc = mm0;
         kloc = 0;
         break;
      case 2:
         iloc = 0;
         jloc = 0;
         kloc = mm0;
         break;
      }
   }

};

}  // namespace BlkC

#endif
