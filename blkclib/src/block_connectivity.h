#ifndef _BLOCK_CONNECTIVITY_INCLUDED
#define _BLOCK_CONNECTIVITY_INCLUDED

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
#include <stdexcept>

#ifdef BLKC_ENABLE_DEBUGGING
#include <fstream>
#endif


/*******************************************************************************
 *
 * Special types
 *
 ******************************************************************************/

namespace BlkC
{

enum ILoc_t {
   IMin = -1,
   IAll = 0,
   IMax = 1
};

enum JLoc_t {
   JMin = -1,
   JAll = 0,
   JMax = 1
};

enum KLoc_t {
   KMin = -1,
   KAll = 0,
   KMax = 1
};

enum Cell_t {
   CellCenter,
   CellVertex
};


/*******************************************************************************
 *
 * class VertexPack
 *
 * Purpose
 * =======
 *
 *   Collects the corner data required of a block and is used to transfer the
 *   information to the BlockConnectivity class.  Multiple methods of loading
 *   data are provided.  The coordinates of each physical (as opposed to ghost)
 *   corner vertex is required.  Vectors, pointing away from the corner, are
 *   also required along each edge that connects to a corner.  E.g., the vector
 *   in the i direction at vertex (IMin, JMin, KMin) is given by
 *   Node(1, 0, 0) - Node(0, 0, 0).  The vectors should describe as closely as
 *   possible the direction of the edge near the corner (going away from the
 *   corner).  The magnitude of the vector is not important.
 *
 *   A vertex can be described by using the enumerations (ILoc_t, JLoc_T,
 *   KLoc_t), a similar set of integers with values -1 or 1, or an
 *   integer m.  See the notes for class BlockConnectivity for how m relates to
 *   the computational indices.
 *
 * Member Functions
 * ================
 *
 *   void set_coord(const int iloc, const int jloc, const int kloc,
 *             const double x, const double y, const double z)
 *                     -- sets the physical coordinate of a corner vertex
 *     iloc, jloc, kloc   (I) corner of the block
 *     x, y, z            (I) physical coordinates of the corner
 *
 *   void set_coord(const int m, const double x, const double y, const double z)
 *                     -- sets the physical coordinate of a corner vertex
 *     m                  (I) Array index of the corner vertex
 *     x, y, z            (I) physical coordinates of the corner
 *
 *   void set_coord(const int iloc, const int jloc, const int kloc,
 *             const double *const crd)
 *                     -- sets the physical coordinate of a corner vertex
 *     iloc, jloc, kloc   (I) corner of the block
 *     crd                (I) rank 1 array of dimension [3] describing the
 *                            physical coordinatex x, y, z
 *
 *   void set_coord(const int m, const double *crd)
 *                     -- sets the physical coordinate of a corner vertex
 *     m                  (I) Array index of the corner vertex
 *     crd                (I) rank 1 array of dimension [3] describing the
 *                            physical coordinatex x, y, z
 *
 *   void set_vector(const int iloc, const int jloc, const int kloc, char dir,
 *                   const double vx, const double vy, const double vz)
 *                     -- sets the edge vector at a corner vertex (pointing out
 *                        from the corner along the edge)
 *     iloc, jloc, kloc   (I) corner of the block
 *     dir                (I) computation direction of the edge (sign not
 *                            important).  Permitted values are 'i', 'j', or
 *                            'k'.
 *     vx, vy, vz         (I) components of the vector
 *
 *   void set_vector(const int m, char dir, const double vx, const double vy,
 *                   const double vz)
 *                     -- sets the edge vector at a corner vertex (pointing out
 *                        from the corner along the edge)
 *     m                  (I) Array index of the corner vertex
 *     dir                (I) computation direction of the edge (sign not
 *                            important).  Permitted values are 'i', 'j', or
 *                            'k'.
 *     vx, vy, vz         (I) components of the vector
 *
 *   void set_vector(const int iloc, const int jloc, const int kloc, char dir,
 *                   const double *const vec)
 *                     -- sets the edge vector at a corner vertex (pointing out
 *                        from the corner along the edge)
 *     iloc, jloc, kloc   (I) corner of the block
 *     dir                (I) computation direction of the edge (sign not
 *                            important).  Permitted values are 'i', 'j', or
 *                            'k'.
 *     vec                (I) rank 1 array of dimension [3] describing the
 *                            components of the vector
 *
 *   void set_vector(const int m, char dir, const double vx, const double vy,
 *                   const double *vec)
 *                     -- sets the edge vector at a corner vertex (pointing out
 *                        from the corner along the edge)
 *     m                  (I) Array index of the corner vertex
 *     dir                (I) computation direction of the edge (sign not
 *                            important).  Permitted values are 'i', 'j', or
 *                            'k'.
 *     vec                (I) rank 1 array of dimension [3] describing the
 *                            components of the vector
 *
 *   int check_flag()  -- checks if all values have been set
 *     return             0   - all values have been set
 *                        1-8 - first vertex found with a value not set
 *
 * Notes
 * =====
 *
 *   - when BlockConnectivity reads a VertexPack, it first checks that data
 *     has been entered for each vertex and then clears the flag.
 *
 ******************************************************************************/


class VertexPack
{

  public:
   VertexPack()
   {
      clear_flag();
   }

   void set_coord(const int iloc, const int jloc, const int kloc,
                  const double x, const double y, const double z)
   {
      set_coord(ijktom(iloc, jloc, kloc), x, y, z);
   }

   void set_coord(const int m, const double x, const double y, const double z)
   {
      double *const p = coord[m];
      p[0] = x;
      p[1] = y;
      p[2] = z;
      flag[m][0] = false;
   }

   void set_coord(const int iloc, const int jloc, const int kloc,
                  const double *const crd)
   {
      set_coord(ijktom(iloc, jloc, kloc), crd);
   }

   void set_coord(const int m, const double *crd)
   {
      double *const p = coord[m];
      p[0] = crd[0];
      p[1] = crd[1];
      p[2] = crd[2];
      flag[m][0] = false;
   }

   void set_vector(const int iloc, const int jloc, const int kloc, char dir,
                   const double vx, const double vy, const double vz)
   {
      set_vector(ijktom(iloc, jloc, kloc), dir, vx, vy, vz);
   }

   void set_vector(const int m, char dir, const double vx, const double vy,
                   const double vz)
   {
      const int c = dirtoc(dir);
      double *const p = vector[m][c];
      const double mag = std::sqrt(vx*vx + vy*vy + vz*vz);
      p[0] = vx/mag;
      p[1] = vy/mag;
      p[2] = vz/mag;
      flag[m][c+1] = false;
   }

   void set_vector(const int iloc, const int jloc, const int kloc, char dir,
                   const double *const vec)
   {
      set_vector(ijktom(iloc, jloc, kloc), dir, vec);
   }

   void set_vector(const int m, char dir, const double *vec)
   {
      const int c = dirtoc(dir);
      double *p = vector[m][c];
      const double mag =
         std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
      p[0] = vec[0]/mag;
      p[1] = vec[1]/mag;
      p[2] = vec[2]/mag;
      flag[m][c+1] = false;
   }

   int check_flag()
   {
      for ( int m = 0; m != 8; ++m ) {
         if ( flag[m][0] || flag[m][1] || flag[m][2] || flag[m][3] )
            return m+1;
      }
      return 0;
   }

  private:
   double coord[8][3];
   double vector[8][3][3];
   bool flag[8][4];                     // T = data not entered

   int ijktomv(const int iloc, const int jloc, const int kloc) const
   {
      return kloc + 1 >> 1 << 2 | jloc + 1 >> 1 << 1 | iloc + 1 >> 1;
   }

   int dirtoc(const char dir) const
   {
      switch ( dir ) {
      case 'i':
         return 0;
      case 'j':
         return 1;
      case 'k':
         return 2;
      default:
         throw std::invalid_argument("Invalid direction given for vector");
      }
   }

   void clear_flag() 
   {
      for ( int m = 0; m != 8; ++m ) {
         flag[m][0] = true;
         flag[m][1] = true;
         flag[m][2] = true;
         flag[m][3] = true;
      }
   }

   // Convert ijk indices to vertex index "mv" in array
   static int ijktom(const int iloc, const int jloc, const int kloc)
   {
      // Note that + 1 >> 1 converts -1 to 0
      return
         kloc + 1 << 1 |
         jloc + 1 |
         iloc + 1 >> 1;
   }

   // Convert vertex index "mv" in array to ijk indices
   static void mtoijk(const int cloc, int &iloc, int &jloc, int &kloc)
   {
      // Note that (<< 1) - 1 converts 0 to -1
      iloc = ((cloc & 1) << 1) - 1;
      jloc = ((cloc >> 1 & 1) << 1) - 1;
      kloc = (cloc >> 2 << 1) - 1;
   }

   friend class BlockConnectivity;

};


/*******************************************************************************
 *
 * Class BlockConnectivity
 *
 * Purpose
 * =======
 *
 *   A class that computes the connectivity between a set of unstructured blocks
 *   as well as the orientation of each block.  Transformation parameters are
 *   provided to compute the transformation of an index from one block to
 *   another.  The transformations and offsets all reference a "primary" block.
 *   This is taken to be the first block added to this library.
 *   - Transformations transform FROM the primary to "this" block.
 *   - Offsets are FROM the primary to "this" block.
 *
 *   See the CGNS docs:
 *   <http://www.grc.nasa.gov/WWW/cgns/sids/cnct.html#GridConnectivity1to1>
 *   for details on the transformation matrix.  In this class, begin1 is always
 *   (0, 0, 0) and begin2 is defined as the offset.
 *
 * Constructors
 * ============
 *
 *   BlockConnectivity(const int numBlk_in, const Cell_t cellType_in,
 *                     const int numGhostLayer_in, const int precision = 6)
 *                     -- initializes the class
 *     numBlock_in        (I) number of blocks in the domain
 *     cellType_in        (I) BlkC::CellCenter for cell centered grids
 *                            BlkC::CellVertex for cell vertex grids
 *     numGhostLayer_in   (I) number of ghost layers.  Assumed to be the same in
 *                            each direction
 *     precision          (I) Precision for matching vertices.  Gives number of
 *                            significant digits for a match between two
 *                            doubles.  Default is 6.
 *
 * Member Functions
 * ================
 *   
 *   void add_block(const int blockIndex,
 *                  const int blkDimI,
 *                  const int blkDimJ,
 *                  const int blkDimK,
 *                  VertexPack &vp)
 *                     -- adds a block to the class.
 *                        *** All blocks must be added before any of the
 *                            routines which query the block connectivity (e.g.,
 *                            'num_neighbour_block', 'neighbour_block',
 *                            'at_domain_extent') are called.
 *     blockIndex         (I) Index of the block to add.  The index directly
 *                            corresponds to a location in an array and must be
 *                            in the range (0 <= blockIndex < numBlock_in) where
 *                            numBlock_in is the number of blocks specified in
 *                            the constructor
 *     blkDimI, blkDimK, blkDimJ
 *                        (I) Dimensions of the block
 *     vp                 (I) Location and vectors at the corner vertices
 *                        (O) Flag reset to 'not-loaded'
 * 
 *   int num_neighbour_block(const int blockIndex, const int iloc,
 *                           const int jloc, const int kloc) const
 *                     -- number of neighbour blocks across a block boundary
 *                        element
 *     blockIndex         (I) index of local block
 *     iloc, jloc, kloc   (I) index of element across which to find neighbour
 *     return             (O) number of neighbours
 *
 *   int neighbour_block(const int blockIndex, const int iloc, const int jloc,
 *                       const int kloc, const int nbrloc,
 *                        int &blockIndexNbr, int &ilocNbr, int &jlocNbr,
 *                        int &klocNbr, int *const ctm) const
 *                     -- data for a neighbour block across a block boundary
 *                        element
 *     blockIndex         (I) index of local block
 *     iloc, jloc, kloc   (I) index of element across which to find neighbour
 *     nbrloc             (I) Neighbour to find.  Start counting from 0.
 *     blockIndexNbr      (O) Index of the neighbour block
 *     ilocNbr, jlocNbr, klocNbr
 *                        (O) Index of the shared boundary element in the
 *                            neighbour block.
 *     ctm                (O) A rank 1 array of of dimension 6
 *                            [0] to [2] the compact transformation
 *                            [3] to [5] the offset
 *     return             (O) 0 - success
 *                            1 - requested neighbour not found
 * 
 *   bool at_domain_extent(const int blockIndex, const int iloc, const int jloc,
 *                        const int kloc) const
 *                     -- indicate whether a block boundary element (vertex,
 *                        edge, or face) is at an extent of the domain.  This is
 *                        useful for determining if block has no neighbours
 *                        across a given boundary element because of
 *                        unstructured connectivity or because it is at the edge
 *                        of the domain.
 *     blockIndex         (I) index of local block
 *     iloc, jloc, kloc   (I) index of element across which to find neighbour
 *     return             (O) T - at domain extent
 *                            F - internal
 *
 * Notes
 * =====
 *
 *   - Conventions for node numbering in a block
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
 *   - All boundary elements (faces, edges, and vertices) are describe using
 *     the notation (i, j, k) where i, j, and k can have the values -1, 0, or 1.
 *     Alternatively, the enumerations ILoc_t, JLoc_t, and KLoc_t can be used
 *     to specify *Min, *All, and *Max.
 *     * faces are indicated by setting one coordinate to non-zero values.  E.g,
 *       face (0, 4, 6, 2) is given by (IMin, JAll, KAll) or (-1, 0, 0).
 *     * edges are indicated by setting two coordinates to non-zero values.
 *       E.g., edge (0, 1) is given by (IAll, JMin, KMin) or (0, -1, -1).
 *     * vertices are indictaed by setting all coordinates to non-zero values.
 *       E.g., vertex (0) is given by (IMin, JMin, KMin) or (-1, -1, -1).
 *   - All blocks must be added before calling one of the routines that queries
 *     the block connectivity.
 *   - Debugging code is defined by BLKC_ENABLE_DEBUGGING.  This is for
 *     development purposes only
 *
 ******************************************************************************/

class BlockConnectivity
{


/*==============================================================================
 * Constructors and destructors
 *============================================================================*/

  public:

   // The constructor
   BlockConnectivity(const int numBlk_in, const Cell_t cellType_in,
                     const int numGhostLayer_in, const int precision = 6);

   // Destructor
   ~BlockConnectivity();

   // Copy and assignment are not permitted


/*==============================================================================
 * Member functions
 *============================================================================*/

   // Add a block
   void add_block(const int blockIndex,
                  const int blkDimI,
                  const int blkDimJ,
                  const int blkDimK,
                  VertexPack &vp);

   // Number of neighbour blocks across a block boundary element
   int num_neighbour_block(const int blockIndex, const int iloc,
                           const int jloc, const int kloc);

   // Data for a neighbour block across a block boundary element
   int neighbour_block(const int blockIndex, const int iloc, const int jloc,
                       const int kloc, const int nbrloc,
                       int &blockIndexNbr, int &ilocNbr, int &jlocNbr,
                       int &klocNbr, int *const ctm);

   // If boundary element is at extent of the domain
   bool at_domain_extent(const int blockIndex, const int iloc, const int jloc,
                         const int kloc);

#ifdef BLKC_ENABLE_DEBUGGING
   void dbg_save_input(const char *const fileName);
   void dbg_load_input();
  private:
   std::ofstream fout;
   int sizeBlkList;
#endif


/*==============================================================================
 * Data members
 *============================================================================*/

  private:
   const Cell_t cellType;               // Type of grid (cell-centered or cell-
                                        // vertex)
   const int numGhostLayer;             // Number of ghost layers (assumed to
                                        // be the same in each direction)
   int numBlk;                          // Number of blocks in the domain.  This
                                        // is the number added, not the number
                                        // specified by the user in the
                                        // constructor
   int uniqueVert;                      // Number of unique vertices in the
                                        // domain
   bool computeCalled;                  // T - compute_neighbours has been
                                        //     called
                                        // F - compute_neighbours has not been
                                        //     called
                                        // compute_neighbours() may only be
                                        // called once

//--The following data types are stored as void* to avoid including the class
//--definitions

   void *_blkList;                      // List of blocks in the domain
   void *_nbrDataMem;                   // Memory for the neighbour data
                                        // structure
   void *_vertexMap;                    // Map for finding unique vertices


/*==============================================================================
 * Private constructors
 *============================================================================*/

   BlockConnectivity(const BlockConnectivity&);
   BlockConnectivity &operator=(const BlockConnectivity&);


/*==============================================================================
 * Private member functions
 *============================================================================*/

   void compute_neighbours();

};

}  // namespace BlkC

#endif
