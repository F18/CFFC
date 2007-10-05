
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

#include <iterator>
#include <map>
#include <set>
#include <stdexcept>

#ifdef HAVE_HASH_MAP
#include "hash_map.h"
#endif

#include "block_connectivity.h"
#include "BlkC_block.h"
#include "BlkC_parameters.h"
#include "BlkC_vertex.h"
#include "mesh_neighbour.h"
#include "Pool_create.h"
#include "Sort.h"
#include "coord_transform.h"

#ifdef BLKC_ENABLE_DEBUGGING
#include <fstream>
#endif


/*==============================================================================
 * Types
 *============================================================================*/

struct VertexRep
{
   int index;
   VertexRep() : index(-1) { }
};

#ifdef HAVE_HASH_MAP
typedef HASH_MAP<BlkC::Vertex, VertexRep, BlkC::Hash_Vertex, BlkC::Equal_Vertex>
VertexMap;
#else
typedef std::map<BlkC::Vertex, VertexRep, BlkC::Less_Vertex> VertexMap;
#endif

typedef Pool::Create<BlkC::NeighbourData> PoolNbrData;

typedef std::set<BlkC::Face, BlkC::Less_Face> FaceSet;


/*******************************************************************************
 *
 * Definition for static const in BlkC::Block
 *
 ******************************************************************************/

// Definitions of edges. Also defines vertex rotation patterns for matching
// coordinated systems.
const int BlkC::Block::mEdge[12][2] = {
   // I Dir. {JMin, KMin}, {JMax, KMin}, {JMin, KMax}, {JMax, KMax}
   {0, 1}, {2, 3}, {4, 5}, {6, 7},
   // J Dir. {IMin, KMin}, {IMax, KMin}, {IMin, KMax}, {IMax, KMax}
   {0, 2}, {1, 3}, {4, 6}, {5, 7},
   // K Dir. {IMin, JMin}, {IMax, JMin}, {IMin, JMax}, {IMax, JMax}
   {0, 4}, {1, 5}, {2, 6}, {3, 7}
};

// Definitions of faces. Also defines vertex rotation patterns for matching
// coordinated systems.
const int BlkC::Block::mFace[6][4] = {
   // I Dir. {IMin}, {IMax}
   {0, 4, 6, 2}, {1, 3, 7, 5},
   // J Dir. {JMin}, {JMax}
   {0, 1, 5, 4}, {2, 6, 7, 3},
   // K Dir. {KMin}, {KMax}
   {0, 2, 3, 1}, {4, 5, 7, 6}
};


/*******************************************************************************
 *
 * Definitions for external variables in "BlkC_parameters.h"
 *
 ******************************************************************************/

double BlkC::toler = 1.E-6;
double BlkC::halftoler = 5.E-7;


/*******************************************************************************
 *
 * Helper routines
 *
 ******************************************************************************/


/*==============================================================================
 *
 * Routine erase_common
 *
 * Purpose
 * =======
 *
 *   Removes the common elements between sorted containers a1 and a2 from
 *   container a1.
 *
 * I/O
 * ===
 *
 *   a                  - (I) Iterator to the start of container CA
 *   aEnd               - (I) End iterator to container CA
 *   b                  - (I) Const iterator to the start of container CB
 *   n1                 - (I) End iterator to container CB
 *   return             - new size of container a
 *
 * Notes
 * =====
 *
 *   - Elements are shifted rather than deleted
 *   - Duplicate common elements in a1 are deleted
 *   - Containers should be sorted by <
 *
 *============================================================================*/

template <typename CA, typename CB>
int erase_common(CA a, const CA aEnd, CB b, const CB bEnd)
{
   typedef typename std::iterator_traits<CA>::value_type T;
   int n = 0;
   CA ap = a;
   while ( a != aEnd && b != bEnd ) {
      const T elem_a = *a;
      const T elem_b = *b;
      if ( elem_a < elem_b ) {
         *ap++ = *a++;
         ++n;
      }
      else if ( elem_b < elem_a )
         ++b;
      else
         ++a;
   }
   while ( a != aEnd ) {
      *ap++ = *a++;
      ++n;
   }
   return n;
}


/*==============================================================================
 *
 * Routine rotate_value
 *
 * Purpose
 * =======
 *
 *   Breaks down a "rotate value" into a sign and an index
 *
 *============================================================================*/

inline void rotate_value(const int rot, int &sign, int &index)
{
   sign = ( rot < 0 ) ? -1 : 1;
   index = std::abs(rot) - 1;
}


/*==============================================================================
 *
 * Routine set_dot
 *
 * Purpose
 * =======
 *
 *   Takes the dot product of two sets of vectors and returns the sum over the
 *   set.  Set 2 can be rotated right by iRot.
 *
 *============================================================================*/

inline double set_dot(const double set1[3][3],
                      const double set2[3][3],
                      const int iRot)
{
   const int s2vi = (3 - iRot) % 3;
   const int s2vj = (4 - iRot) % 3;
   const int s2vk = (5 - iRot) % 3;
   return
      set1[0][0]*set2[s2vi][0] + set1[0][1]*set2[s2vi][1] +
      set1[0][2]*set2[s2vi][2] +
      set1[1][0]*set2[s2vj][0] + set1[1][1]*set2[s2vj][1] +
      set1[1][2]*set2[s2vj][2] +
      set1[2][0]*set2[s2vk][0] + set1[2][1]*set2[s2vk][1] +
      set1[2][2]*set2[s2vk][2];
}


/*******************************************************************************
 *
 * Member routines to class BlockConnectivity
 *
 ******************************************************************************/

/*==============================================================================
 *
 * Routine BlockConnectivity
 *
 * Purpose
 * =======
 *
 *   Constructor
 *
 *============================================================================*/

BlkC::BlockConnectivity::BlockConnectivity(const int numBlk_in,
                                           const Cell_t cellType_in,
                                           const int numGhostLayer_in,
                                           const int precision)
   :
   cellType(cellType_in), numGhostLayer(numGhostLayer_in), numBlk(0),
   uniqueVert(0), computeCalled(false)
{

   _blkList = new Block[numBlk_in];
   _vertexMap = new VertexMap();
   _nbrDataMem = new PoolNbrData(100);
   toler = std::pow(10., -precision);
   halftoler = 0.5*toler;

#ifdef BLKC_ENABLE_DEBUGGING
   sizeBlkList = numBlk_in;
#endif

}


/*==============================================================================
 *
 * Routine ~BlockConnectivity
 *
 * Purpose
 * =======
 *
 *   Destructor
 *
 *============================================================================*/

BlkC::BlockConnectivity::~BlockConnectivity()
{

   Block *blkList = static_cast<Block*>(_blkList);
   delete[] blkList;
   PoolNbrData *nbrDataMem = static_cast<PoolNbrData*>(_nbrDataMem);
   delete nbrDataMem;
   VertexMap *vertexMap = static_cast<VertexMap*>(_vertexMap);
   delete vertexMap;

#ifdef BLKC_ENABLE_DEBUGGING
   if ( fout.is_open() ) fout.close();
#endif

}


/*==============================================================================
 *
 * Routine add_block
 *
 * Purpose
 * =======
 *
 *   Adds a block to the library
 *
 * Notes
 * =====
 *
 *   - See notes in "block_connectivity.h" for node numbering conventions
 *
 *============================================================================*/

void BlkC::BlockConnectivity::add_block(const int blockIndex,
                                        const int blkDimI,
                                        const int blkDimJ,
                                        const int blkDimK,
                                        VertexPack &vp)
{

   if ( computeCalled )
      throw std::runtime_error("Blocks may not be added to class "
                               "BlockConnectivity after calling a routine that"
                               "queries the connectivity.");

//--Required pointers

   Block *const blkList = static_cast<Block*>(_blkList);
   VertexMap *const vertexMap = static_cast<VertexMap*>(_vertexMap);

//--Check the input vp

   if ( vp.check_flag() )
      throw std::invalid_argument("All required data has not been set in "
                                  "class VertexPackage.");
   vp.clear_flag();

//--Describe the vertices by unique integers rather than spatial coordinates in
//--the block list

   Vertex vert_coord[8];
   vert_coord[0].set_coord(vp.coord[0]);
   vert_coord[1].set_coord(vp.coord[1]);
   vert_coord[2].set_coord(vp.coord[2]);
   vert_coord[3].set_coord(vp.coord[3]);
   vert_coord[4].set_coord(vp.coord[4]);
   vert_coord[5].set_coord(vp.coord[5]);
   vert_coord[6].set_coord(vp.coord[6]);
   vert_coord[7].set_coord(vp.coord[7]);

   Block &block = blkList[blockIndex];
   block.label = blockIndex;
   block.dim[0] = blkDimI - 1;
   block.dim[1] = blkDimJ - 1;
   block.dim[2] = blkDimK - 1;
   for ( int i = 0; i != 8; ++i ) {
      VertexRep &vrep = (*vertexMap)[vert_coord[i]];
      if ( vrep.index < 0 ) vrep.index = uniqueVert++;
      block.set_vertex(i, vrep.index);
   }
   block.set_vector(vp.vector[0][0]);
   ++numBlk;

#ifdef BLKC_ENABLE_DEBUGGING
   if ( fout.is_open() ) {
      fout.setf(std::ios_base::scientific, std::ios_base::floatfield);
      fout.precision(std::numeric_limits<double>::digits10);
      fout << std::endl << blockIndex << std::endl;
      for ( int i = 0; i != 8; ++i ) {
         fout << block.get_vertex(i) << ' ' << vp.coord[i][0] << ' '
              << vp.coord[i][1] << ' ' << vp.coord[i][2] << std::endl;
         for ( int j = 0; j != 3; ++j )
            fout << "  " << vp.vector[i][j][0] << ' ' << vp.vector[i][j][1]
                 << ' ' << vp.vector[i][j][2] << std::endl;
      }
      fout.setf(std::ios_base::fmtflags(0), std::ios_base::floatfield);
      fout.precision(6);
   }
#endif

}


/*==============================================================================
 *
 * Routine compute_neighbours
 *
 * Purpose
 * =======
 *
 *   After all blocks have been added, this routine is called to:
 *   a) create a database of the block connectivity
 *   b) compute the transformation and offset from the prime reference block
 *      (the first block) to each other block
 *   c) store information about each blocks neighbours including the label of
 *      the neighbour block and which boundary is shared with "this" block.
 *
 * Notes
 * =====
 *
 *   - Neighbour blocks are often referred to as donor blocks
 *
 *============================================================================*/

//--Description of boundary (vertex, edge, or face) indices

struct Indices
{
   int m;
   int i;
   int j;
   int k;
};

void BlkC::BlockConnectivity::compute_neighbours()
{

   if ( computeCalled )
      throw std::runtime_error("compute_neighbours() may only be called once.");
   computeCalled = true;

//--Required pointers

   Block *const blkList = static_cast<Block*>(_blkList);
   PoolNbrData *nbrDataMem = static_cast<PoolNbrData*>(_nbrDataMem);

   MNeighbour meshNeighbour;
   FaceSet faceSet;

//--Add the blocks to the database

   for ( int iBlk = 0; iBlk != numBlk; ++iBlk )
      if ( blkList[iBlk].label != -1 )
         meshNeighbour.add_block(&blkList[iBlk]);

//--Find maximum blocks connected to vertices, edges, and faces

   // These are the inclusive maximums (all blocks connected to vertex, edge, or
   // face)
   const int maxVertNbr = meshNeighbour.max_vertex_neighbours();
   const int maxEdgeNbr = meshNeighbour.max_edge_neighbours();
   const int maxFaceNbr = meshNeighbour.max_face_neighbours();

   // Conservative estimates for the number of blocks connected to a given block
   // via vertices, edges, and faces.
   const int maxBlkEdgeNbr = 12*(maxEdgeNbr - 1);
   const int maxBlkFaceNbr = 6*(maxFaceNbr - 1);

   // Allocate space for receiving or storing neighbour blocks
   Indices *nbrIndices = new Indices[maxVertNbr];
                                        // Storage of indices for neighbour
                                        // blocks
   Block **nbrList = new Block*[maxVertNbr];
                                        // List of neighbour blocks across an
                                        // edge or vertex
   Block **nbrFaceList = new Block*[maxBlkFaceNbr];
                                        // List of neighbour blocks that connect
                                        // to the local block uniquely across
                                        // the faces
   Block **nbrEdgeList = new Block*[maxBlkEdgeNbr + maxBlkFaceNbr];
                                        // List of neighbour blocks that connect
                                        // to the local block uniquely across
                                        // the edges

//--Loop over all blocks

   for ( int iBlk = 0; iBlk != numBlk; ++iBlk ) {
      Block *localBlock = &blkList[iBlk];
      if ( localBlock->label != -1 ) {

         // All neighbour blocks that connect via faces, edges, and vertices
         const int numNbrAllFace =
            meshNeighbour.all_block_face_neighbours(localBlock, nbrFaceList);
         Sort::intro(numNbrAllFace, nbrFaceList);
         const int numNbrAllEdge =
            meshNeighbour.all_block_edge_neighbours(localBlock, nbrEdgeList);
         Sort::intro(numNbrAllEdge, nbrEdgeList);


/*==============================================================================
 *
 * Face Neighbours
 * ===============
 *
 * See the notes section for the reference block and coordinate system.
 *
 * The local block is orientated with the local face at imin of the reference
 * block and such that the other coordinates increase in the reference j and k
 * directions.  CTMRefToLocalFace describes exactly how the local coordinates
 * align with the reference coordindates.
 *
 * The donor block is orientated with the donor face at imax of the reference
 * block and such that the other coordinates increase in the reference j and k
 * directions.  CTMRefToDonorFace describes exactly how the donor coordinates
 * align with the reference coordinates.
 *
 * The donor block can rotate and still match faces with the local block.  The
 * amount of rotation is specified by matching a reference vertex in the local
 * block (mFace[face][0]) with a vertex in the donor block (mFace[face][*]).
 * CTMRefToDonorFace information can still be used but must be rotated through a
 * specific pattern as given by rotCTMFace.  rotCTMFace basically provides
 * indices to CTMRefToDonorFace.  Note that rotCTMFace[0][*] provides the
 * original specification of CTMRefToDonorFace.
 *
 * The offsets are defined in the coordinate system of the reference frame.  The
 * offsets always point from vertex 0 in the block to the reference vertex 0.
 * When defining offsets, the donor block must be placed in its actual position
 * (in the reference frame) with respect to the local block.  Considerations are
 * made for ghost cells and cell-centered meshes.  The required offset points
 * from vertex 0 in the donor block to vertex 0 in the local block.
 *
 *============================================================================*/

         // For each face of the block
         for ( int mLocalFace = 0; mLocalFace != 6; ++mLocalFace ) {

            // Coordinates of the local face described by ijk indices
            int iLocalFace, jLocalFace, kLocalFace;
            Block::mftoijk(mLocalFace, iLocalFace, jLocalFace, kLocalFace);

            // Check if this face neighbour has yet been set (if the neighbour
            // is set, then the extents are also known)
            if ( !localBlock->have_neighbour(iLocalFace, jLocalFace,
                                             kLocalFace) ) {

               // Representation of the face
               Face face = localBlock->get_face(mLocalFace);
               // Find the block that connects across this face
               Block *donorBlock = meshNeighbour.block_face_neighbour(
                  localBlock, face);
               if ( donorBlock ) {  // Neighbour exists

                  // Coordinates of the donor face described by ijk indices
                  const int mDonorFace = donorBlock->find_face_index(face);
                  int iDonorFace, jDonorFace, kDonorFace;
                  Block::mftoijk(mDonorFace,
                                 iDonorFace, jDonorFace, kDonorFace);

                  // Set extents for the face
                  localBlock->set_extent(iLocalFace, jLocalFace, kLocalFace,
                                         false);
                  donorBlock->set_extent(iDonorFace, jDonorFace, kDonorFace,
                                         false);

//--Describe the transformation

                  // Compact transformation matrices describing the orientation
                  // of the local block with respect to the reference frame
                  static const int CTMRefToLocalFace[6][3] = {
                     { 1,  2,  3}, {-1,  3,  2}, { 2,  3,  1},
                     {-2,  1,  3}, { 3,  1,  2}, {-3,  2,  1}
                  };

                  // Compact transformation matrices describing the orientation
                  // of the donor block with respect to the reference frame
                  static const int CTMRefToDonorFace[6][3] = {
                     {-1,  3,  2}, { 1,  2,  3}, {-2,  1,  3},
                     { 2,  3,  1}, {-3,  2,  1}, { 3,  1,  2}
                  };

                  // Adjustments to CTMRefToDonorFace based on rotation of the
                  // donor block
                  static const int rotCTMFace[4][3] = {
                     { 1,  2,  3}, { 1,  3, -2}, { 1, -2, -3}, { 1, -3,  2}
                  };

                  int vRef =            // Reference vertex
                     localBlock->get_vertex(Block::mFace[mLocalFace][0]);
                  int iRot = 0;         // Amount of rotation in donor block
                  int rotIndex;         // Index into CTMRefToDonorFace as
                                        // modified by rotation
                  int rotSign;          // Sign of CTMRefToDonorFace as modified
                                        // by rotation

                  // Find the amount of rotation in the donor block
                  for ( ; iRot != 4; ++iRot )
                     if ( vRef == donorBlock->get_vertex(
                             Block::mFace[mDonorFace][iRot]) ) break;
                  if ( iRot == 4 )
                     throw std::runtime_error("Error matching vertices between "
                                              "faces even though faces "
                                              "supposedly match.");
                  int ctm[6];           // Compact transformation matrix

                  // Transformations between the local and reference
                  ctm[0] = CTMRefToLocalFace[mLocalFace][0];
                  ctm[1] = CTMRefToLocalFace[mLocalFace][1];
                  ctm[2] = CTMRefToLocalFace[mLocalFace][2];
                  ctm[3] = 0;
                  ctm[4] = 0;
                  ctm[5] = 0;

                  CoordTransform TMRefToLocal(ctm);
                  CoordTransform TMLocalToRef(TMRefToLocal);
                  TMLocalToRef.transpose();

                  // Transformations between the donor and reference
                  rotate_value(rotCTMFace[iRot][0], rotSign, rotIndex);
                  ctm[0] = rotSign*CTMRefToDonorFace[mDonorFace][rotIndex];
                  rotate_value(rotCTMFace[iRot][1], rotSign, rotIndex);
                  ctm[1] = rotSign*CTMRefToDonorFace[mDonorFace][rotIndex];
                  rotate_value(rotCTMFace[iRot][2], rotSign, rotIndex);
                  ctm[2] = rotSign*CTMRefToDonorFace[mDonorFace][rotIndex];

                  CoordTransform TMRefToDonor(ctm);
                  CoordTransform TMLocalToDonor(TMRefToDonor, TMLocalToRef);

//--Describe the offset

                  // Offset in an index caused by rotation of the donor block
                  static const int rotOffsetFace[4][3] = {
                     { 0,  0,  0}, { 0,  0, -1}, { 0, -1, -1}, { 0, -1,  0}
                  };

                  // Offset for local block.  ijk maximum faces invoke an
                  // increment in the -i direcion.  I.e., localOffsetFace[6][3]
                  // is:
                  // { { 0,  0,  0}, {-1,  0,  0}, { 0,  0,  0},
                  //   {-1,  0,  0}, { 0,  0,  0}, {-1,  0,  0} }
                  int localOffset[3];
                  localOffset[0] =
                     -(std::max(0, iLocalFace) +
                       std::max(0, jLocalFace) +
                       std::max(0, kLocalFace));
                  localOffset[1] = 0;
                  localOffset[2] = 0;

                  // The localOffset, defined in the reference frame, needs to
                  // be transformed into the local frame, just so it can be
                  // multiplied by the block dimensions.  Then it needs to be
                  // tranformed to the donor frame.
                  TMRefToLocal.transform(localOffset);  // To local frame
                  localOffset[0] *= localBlock->dim[0];
                  localOffset[1] *= localBlock->dim[1];
                  localOffset[2] *= localBlock->dim[2];
                  // Considerations for ghost cells
                  localOffset[0] += numGhostLayer;
                  localOffset[1] += numGhostLayer;
                  localOffset[2] += numGhostLayer;
                  TMLocalToDonor.transform(localOffset);  // Put in donor frame

                  // Offset for donor block.  ijk maximum faces invoke an
                  // increment in the +i direction.  I.e., donorOffsetFace[6][3]
                  // is:
                  // { { 0,  0,  0}, { 1,  0,  0}, { 0,  0,  0},
                  //   { 1,  0,  0}, { 0,  0,  0}, { 1,  0,  0} }
                  int donorOffset[3];
                  donorOffset[0] =
                     std::max(0, iDonorFace) +
                     std::max(0, jDonorFace) +
                     std::max(0, kDonorFace);

                  // Offset caused by rotation of the donor block
                  // donorOffset[0] += rotOffsetFace[iRot][0];  (Always 0)
                  donorOffset[1] = rotOffsetFace[iRot][1];
                  donorOffset[2] = rotOffsetFace[iRot][2];

                  // The donorOffset, defined in the reference frame, needs to
                  // be transformed into the donor frame and multiplied by the
                  // block dimensions.
                  TMRefToDonor.transform(donorOffset);  // To donor frame
                  donorOffset[0] *= donorBlock->dim[0];
                  donorOffset[1] *= donorBlock->dim[1];
                  donorOffset[2] *= donorBlock->dim[2];
                  // Considerations for ghost cells
                  donorOffset[0] += numGhostLayer;
                  donorOffset[1] += numGhostLayer;
                  donorOffset[2] += numGhostLayer;

                  // For cell-center transformations, add a jump across the
                  // boundary
                  if ( cellType == CellCenter ) {
                     donorOffset[0] += iDonorFace;
                     donorOffset[1] += jDonorFace;
                     donorOffset[2] += kDonorFace;
                  }

                  // Subtract the local_to_ref offset from the donor_to_ref
                  // offset to get the total donor_to_local offset
                  int totalOffset[3];
                  totalOffset[0] = donorOffset[0] - localOffset[0];
                  totalOffset[1] = donorOffset[1] - localOffset[1];
                  totalOffset[2] = donorOffset[2] - localOffset[2];
                  TMLocalToDonor.set_offset(totalOffset);

                  // Get the reverse transformation
                  CoordTransform TMDonorToLocal(TMLocalToDonor);
                  TMDonorToLocal.reverse();

//--Set the neighbour information

                  // Set the donor information as seen from the local block
                  NeighbourData *nbrData = nbrDataMem->get();
                  nbrData->label = donorBlock->label;
                  nbrData->boundary[0] = iDonorFace;
                  nbrData->boundary[1] = jDonorFace;
                  nbrData->boundary[2] = kDonorFace;
                  TMLocalToDonor.compact(nbrData->ctm);
                  localBlock->set_neighbour_data(iLocalFace, jLocalFace,
                                                 kLocalFace, nbrData);

                  // Set the local information as seen from the donor block
                  nbrData = nbrDataMem->get();
                  nbrData->label = localBlock->label;
                  nbrData->boundary[0] = iLocalFace;
                  nbrData->boundary[1] = jLocalFace;
                  nbrData->boundary[2] = kLocalFace;
                  TMDonorToLocal.compact(nbrData->ctm);
                  donorBlock->set_neighbour_data(iDonorFace, jDonorFace,
                                                 kDonorFace, nbrData);

               }  // End if donor block exists
               else {
                  localBlock->set_extent(iLocalFace, jLocalFace, kLocalFace,
                                         true);
               }
            }  // End if this neighbour data not yet set
         }  // End of loop over local faces


/*==============================================================================
 *
 * Edge neighbours
 * ===============
 *
 * See the notes section for the reference block and coordinate system.
 *
 * The local block is orientated with the local edge at edge (iall,jmin,kmin) of
 * the reference block and such that the coordinate increases in the reference
 * i direction.  CTMRefToLocalEdge describes exactly how the local coordinates
 * align with the reference coordinates.
 *
 * The donor block is orientated with the donor edge at edge (iall,jmax,kmax) of
 * the reference block and such that the coordinate increases in the reference i
 * direction.  CTMRefToDonorEdge describes exactly how the donor coordinates
 * align with the reference coordinates.
 *
 * The donor block can swap directions and still match edges with the local
 * block.  The swap is specified by matching a reference vertex in the local
 * block (mEdge[edge][0]) with a vertex in the donor block (mEdge[edge][0|1].
 * CTMRefToDonorEdge information can still be used but must be rotated (swapped)
 * through the specific pattern as given by rotCTMEdge.  rotCTMEdge basically
 * provides indices to CTMRefToDonorEdge.  Note that rotCTMEdge[0][*] provides
 * the original specification of CTMRefToDonorEdge.
 *
 * The offsets are defined in the coordinate system of the reference frame.  The
 * offsets always point from vertex 0 in the block to the reference vertex 0.
 * When defining offsets, the donor block must be placed in its actual position
 * (in the reference frame) with respect to the local block.  Considerations are
 * made for ghost cells and cell-centered meshes.  The required offset points
 * from vertex 0 in the donor block to vertex 0 in the local block.
 *
 *============================================================================*/

         // For each edge of the block
         for ( int mLocalEdge = 0; mLocalEdge != 12; ++mLocalEdge ) {

            // Coordinates of the local edge described by ijk indices
            int iLocalEdge, jLocalEdge, kLocalEdge;
            Block::metoijk(mLocalEdge, iLocalEdge, jLocalEdge, kLocalEdge);

            // Representation of the edge
            Edge edge = localBlock->get_edge(mLocalEdge);
            // Find the blocks that connect across this edge
            int numNbrEdge = meshNeighbour.block_edge_neighbours(
               localBlock, edge, nbrList);

//--Check if at domain extent if not yet known

            if ( !localBlock->extent_known(iLocalEdge, jLocalEdge,
                                           kLocalEdge) ) {
               // Faces from local block
               if ( iLocalEdge ) faceSet.insert(
                  localBlock->get_face(0, jLocalEdge, kLocalEdge));
               if ( jLocalEdge ) faceSet.insert(
                  localBlock->get_face(iLocalEdge, 0, kLocalEdge));
               if ( kLocalEdge ) faceSet.insert(
                  localBlock->get_face(iLocalEdge, jLocalEdge, 0));
               // Faces from neighbour blocks
               for ( int iNbrBlk = 0; iNbrBlk != numNbrEdge; ++iNbrBlk ) {
                  // Coordinates of the donor edge described by ijk inidices
                  nbrIndices[iNbrBlk].m =
                     nbrList[iNbrBlk]->find_edge_index(edge);
                  Block::metoijk(nbrIndices[iNbrBlk].m, nbrIndices[iNbrBlk].i,
                                 nbrIndices[iNbrBlk].j, nbrIndices[iNbrBlk].k);
                  if ( nbrIndices[iNbrBlk].i ) faceSet.insert(
                     nbrList[iNbrBlk]->get_face(
                        0, nbrIndices[iNbrBlk].j, nbrIndices[iNbrBlk].k));
                  if ( nbrIndices[iNbrBlk].j ) faceSet.insert(
                     nbrList[iNbrBlk]->get_face(
                        nbrIndices[iNbrBlk].i, 0, nbrIndices[iNbrBlk].k));
                  if ( nbrIndices[iNbrBlk].k ) faceSet.insert(
                     nbrList[iNbrBlk]->get_face(
                        nbrIndices[iNbrBlk].i, nbrIndices[iNbrBlk].j, 0));
               }
               bool atDomainExtent = false;
               const FaceSet::const_iterator itFaceEnd = faceSet.end();
               for ( FaceSet::const_iterator itFace = faceSet.begin();
                     itFace != itFaceEnd; ++itFace ) {
                  if ( meshNeighbour.face_num_neighbours(*itFace) != 2 ) {
                     atDomainExtent = true;
                     break;
                  }
               }
               faceSet.clear();
               localBlock->set_extent(iLocalEdge, jLocalEdge, kLocalEdge,
                                      atDomainExtent);
               for ( int iNbrBlk = 0; iNbrBlk != numNbrEdge; ++iNbrBlk )
                  nbrList[iNbrBlk]->set_extent(nbrIndices[iNbrBlk].i,
                                               nbrIndices[iNbrBlk].j,
                                               nbrIndices[iNbrBlk].k,
                                               atDomainExtent);
            }

//--Only compute transformations for blocks unique to the edge
               
            Sort::intro(numNbrEdge, nbrList);
            numNbrEdge = erase_common(nbrList, nbrList + numNbrEdge,
                                      nbrFaceList, nbrFaceList + numNbrAllFace);

            for ( int iNbrBlk = 0; iNbrBlk != numNbrEdge; ++iNbrBlk ) {
               Block *donorBlock = nbrList[iNbrBlk];

               // Check if this edge neighbour has yet been set
               if ( !localBlock->have_neighbour(
                       iLocalEdge, jLocalEdge, kLocalEdge,
                       donorBlock->label) ) {

                  // Coordinates of the donor edge described by ijk inidices
                  const int mDonorEdge = donorBlock->find_edge_index(edge);
                  int iDonorEdge, jDonorEdge, kDonorEdge;
                  Block::metoijk(mDonorEdge,
                                 iDonorEdge, jDonorEdge, kDonorEdge);

//--Describe the transformation

                  // Compact transformation matrices describing the orientation
                  // of the local block with respect to the reference frame
                  static const int CTMRefToLocalEdge[12][3] = {
                     { 1,  2,  3}, { 1,  3, -2}, { 1, -3,  2}, { 1, -2, -3},
                     { 2,  3,  1}, { 2, -1,  3}, { 2,  1, -3}, { 2, -3, -1},
                     { 3,  1,  2}, { 3,  2, -1}, { 3, -2,  1}, { 3, -1, -2}
                  };

                  // Compact transformation matrices describing the orientation
                  // of the donor block with respect to the reference frame
                  static const int CTMRefToDonorEdge[12][3] = {
                     { 1, -2, -3}, { 1, -3,  2}, { 1,  3, -2}, { 1,  2,  3},
                     { 2, -3, -1}, { 2,  1, -3}, { 2, -1,  3}, { 2,  3,  1},
                     { 3, -1, -2}, { 3, -2,  1}, { 3,  2, -1}, { 3,  1,  2}
                  };

                  // Adjustments to CTMRefToDonorEdge based on rotation of the
                  // donor block
                  static const int rotCTMEdge[2][3] = {
                     { 1,  2,  3}, {-1,  3,  2}
                  };

                  int vRef =            // Reference vertex
                     localBlock->get_vertex(Block::mEdge[mLocalEdge][0]);
                  int iRot;             // Amount of rotation in the donor block
                  int rotIndex;         // Index into CTMRefToDonorEdge as
                                        // modified by rotation
                  int rotSign;          // Sign of CTMRefToDonorEdge as modified
                                        // by rotation

               // Find the amount of rotation in the donor block
                  if (
                     vRef ==
                     donorBlock->get_vertex(Block::mEdge[mDonorEdge][0]) )
                     iRot = 0;
                  else if (
                     vRef ==
                     donorBlock->get_vertex(Block::mEdge[mDonorEdge][1]) )
                     iRot = 1;
                  else
                     throw std::runtime_error("Error matching vertices between "
                                              "edges even though edges "
                                              "supposedly match.");
                  int ctm[6];           // Compact transformation matrix

                  // Transformations between the local and reference
                  ctm[0] = CTMRefToLocalEdge[mLocalEdge][0];
                  ctm[1] = CTMRefToLocalEdge[mLocalEdge][1];
                  ctm[2] = CTMRefToLocalEdge[mLocalEdge][2];
                  ctm[3] = 0;
                  ctm[4] = 0;
                  ctm[5] = 0;

                  CoordTransform TMRefToLocal(ctm);
                  CoordTransform TMLocalToRef(TMRefToLocal);
                  TMLocalToRef.transpose();

                  // Transformations between the donor and reference
                  rotate_value(rotCTMEdge[iRot][0], rotSign, rotIndex);
                  ctm[0] = rotSign*CTMRefToDonorEdge[mDonorEdge][rotIndex];
                  rotate_value(rotCTMEdge[iRot][1], rotSign, rotIndex);
                  ctm[1] = rotSign*CTMRefToDonorEdge[mDonorEdge][rotIndex];
                  rotate_value(rotCTMEdge[iRot][2], rotSign, rotIndex);
                  ctm[2] = rotSign*CTMRefToDonorEdge[mDonorEdge][rotIndex];

                  CoordTransform TMRefToDonor(ctm);
                  CoordTransform TMLocalToDonor(TMRefToDonor, TMLocalToRef);

//--Describe the offset.

                  // Offset from local block vertex 0 to reference frame 0.
                  // Note:  Offset from donor block vertex 0 to reference frame
                  // 0 is given by negating these values
                  static const int OffsetLocalToRefEdge[12][3] = {
                     { 0,  0,  0}, { 0,  0, -1}, {0, -1,  0}, { 0, -1, -1},
                     { 0,  0,  0}, { 0, -1,  0}, {0,  0, -1}, { 0, -1, -1},
                     { 0,  0,  0}, { 0,  0, -1}, {0, -1,  0}, { 0, -1, -1}
                  };

                  // Offset for local block.
                  int localOffset[3];
                  localOffset[0] = OffsetLocalToRefEdge[mLocalEdge][0];
                  localOffset[1] = OffsetLocalToRefEdge[mLocalEdge][1];
                  localOffset[2] = OffsetLocalToRefEdge[mLocalEdge][2];

                  // The localOffset, defined in the reference frame, needs to
                  // be transformed into the local frame, just so it can be
                  // multiplied by the block dimensions.  Then it needs to be
                  // tranformed to the donor frame.
                  TMRefToLocal.transform(localOffset);  // To local frame
                  localOffset[0] *= localBlock->dim[0];
                  localOffset[1] *= localBlock->dim[1];
                  localOffset[2] *= localBlock->dim[2];
                  // Considerations for ghost cells
                  localOffset[0] += numGhostLayer;
                  localOffset[1] += numGhostLayer;
                  localOffset[2] += numGhostLayer;
                  TMLocalToDonor.transform(localOffset);  // Put in donor frame

                  // Offset for donor block.  This is the same as the negation
                  // of OffsetLocalToRefEdge.  If swapped (iRot = 1) then the
                  // pattern is from (0, b, c) to (-1, c, b).
                  int donorOffset[3];
                  switch ( iRot ) {
                  case 0:
                     donorOffset[0] = 0;
                     donorOffset[1] = -OffsetLocalToRefEdge[mDonorEdge][1];
                     donorOffset[2] = -OffsetLocalToRefEdge[mDonorEdge][2];
                     break;
                  case 1:
                     donorOffset[0] = -1;
                     donorOffset[1] = -OffsetLocalToRefEdge[mDonorEdge][2];
                     donorOffset[2] = -OffsetLocalToRefEdge[mDonorEdge][1];
                     break;
                  }

                  // The donorOffset, defined in the reference frame, needs to
                  // be transformed into the donor frame and multiplied by the
                  // block dimensions.
                  TMRefToDonor.transform(donorOffset);  // To donor frame
                  donorOffset[0] *= donorBlock->dim[0];
                  donorOffset[1] *= donorBlock->dim[1];
                  donorOffset[2] *= donorBlock->dim[2];
                  // Considerations for ghost cells
                  donorOffset[0] += numGhostLayer;
                  donorOffset[1] += numGhostLayer;
                  donorOffset[2] += numGhostLayer;

                  // For cell-center transformations, add a jump across the
                  // boundary
                  if ( cellType == CellCenter ) {
                     donorOffset[0] += iDonorEdge;
                     donorOffset[1] += jDonorEdge;
                     donorOffset[2] += kDonorEdge;
                  }

                  // Subtract the local_to_ref offset from the donor_to_ref
                  // offset to get the total donor_to_local offset
                  int totalOffset[3];
                  totalOffset[0] = donorOffset[0] - localOffset[0];
                  totalOffset[1] = donorOffset[1] - localOffset[1];
                  totalOffset[2] = donorOffset[2] - localOffset[2];
                  TMLocalToDonor.set_offset(totalOffset);

                  // Get the reverse transformation
                  CoordTransform TMDonorToLocal(TMLocalToDonor);
                  TMDonorToLocal.reverse();

//--Set the neighbour information

                  // Set the donor information as seen from the local block
                  NeighbourData *nbrData = nbrDataMem->get();
                  nbrData->label = donorBlock->label;
                  nbrData->boundary[0] = iDonorEdge;
                  nbrData->boundary[1] = jDonorEdge;
                  nbrData->boundary[2] = kDonorEdge;
                  TMLocalToDonor.compact(nbrData->ctm);
                  localBlock->set_neighbour_data(iLocalEdge, jLocalEdge,
                                                 kLocalEdge, nbrData);

                  // Set the local information as seen from the donor block
                  nbrData = nbrDataMem->get();
                  nbrData->label = localBlock->label;
                  nbrData->boundary[0] = iLocalEdge;
                  nbrData->boundary[1] = jLocalEdge;
                  nbrData->boundary[2] = kLocalEdge;
                  TMDonorToLocal.compact(nbrData->ctm);
                  donorBlock->set_neighbour_data(iDonorEdge, jDonorEdge,
                                                 kDonorEdge, nbrData);

               }  // End if this neighbour data not yet set
            }  // Loop over neighbour blocks across this edge
         }  // End of loop over local edges


/*==============================================================================
 *
 * Vertex neighbours
 * =================
 *
 * The local block is orientated with the local vertex at the (imin,jmin,kmin)
 * vertex of the reference block and such that the i coordinate of the local
 * block is aligned with the i coordinate of the reference block.
 *
 * The donor block is orientated with the donor vertex at the (imax,jmax,kmax)
 * vertex of the reference block and such that the i coordinate of the donor
 * block is aligned with the i coordinate of the reference block.
 *
 * The donor block can rotate to 3 orientations about the shared vertex.  The
 * orientation is determined by comparing the aligned by of the edge vectors
 * near the shared vertex between the local and donor blocks.  If the an
 * edge matches exactly, the dot product between the unit vectors that describe
 * the edge direction (for the local and donor blocks) should be -1.  The sum
 * over all three edges gives a weighting and the minimum weighting over all 3
 * orientations is the best orientation.  rotCTMVert provides indices to
 * CTMRefToDonorVert (Note that CTMRefToDonorVert is also given by
 * CTMRefToLocalVert; see below).  Note that rotCTMVert[0][*} provides the
 * original specification of CTMRefToDonorVert.
 *
 * The offsets are defined in the coordinate system of the reference frame.  The
 * offsets always point from vertex 0 in the block to the reference vertex 0.
 * When defining offsets, the donor block must be placed in its actual position
 * (in the reference frame) with respect to the local block.  Considerations are
 * made for ghost cells and cell-centered meshes.  The required offset points
 * from vertex 0 in the donor block to vertex 0 in the local block.
 *
 *============================================================================*/

         // For each vertex of the block
         for ( int mLocalVert = 0; mLocalVert != 8; ++mLocalVert ) {

            // Coordinates of the local vertex described by ijk indices
            int iLocalVert, jLocalVert, kLocalVert;
            Block::mvtoijk(mLocalVert, iLocalVert, jLocalVert, kLocalVert);

            // Representation of the vertex
            int vert = localBlock->get_vertex(mLocalVert);
            // Find the blocks that connect across this vertex
            int numNbrVert = meshNeighbour.block_vertex_neighbours(
               localBlock, vert, nbrList);

//--Check if at domain extent if not yet known

            if ( !localBlock->extent_known(iLocalVert, jLocalVert,
                                           kLocalVert) ) {
               // Faces from local block
               faceSet.insert(localBlock->get_face(iLocalVert, 0, 0));
               faceSet.insert(localBlock->get_face(0, jLocalVert, 0));
               faceSet.insert(localBlock->get_face(0, 0, kLocalVert));
               // Faces from neighbour blocks
               for ( int iNbrBlk = 0; iNbrBlk != numNbrVert; ++iNbrBlk ) {
                  // Coordinates of the donor vertex described by ijk inidices
                  nbrIndices[iNbrBlk].m =
                     nbrList[iNbrBlk]->find_vertex_index(vert);
                  Block::mvtoijk(nbrIndices[iNbrBlk].m, nbrIndices[iNbrBlk].i,
                                 nbrIndices[iNbrBlk].j, nbrIndices[iNbrBlk].k);
                  faceSet.insert(nbrList[iNbrBlk]->get_face(
                                    nbrIndices[iNbrBlk].i, 0, 0));
                  faceSet.insert(nbrList[iNbrBlk]->get_face(
                                    0, nbrIndices[iNbrBlk].j, 0));
                  faceSet.insert(nbrList[iNbrBlk]->get_face(
                                    0, 0, nbrIndices[iNbrBlk].k));
               }
               bool atDomainExtent = false;
               const FaceSet::const_iterator itFaceEnd = faceSet.end();
               for ( FaceSet::const_iterator itFace = faceSet.begin();
                     itFace != itFaceEnd; ++itFace ) {
                  if ( meshNeighbour.face_num_neighbours(*itFace) != 2 ) {
                     atDomainExtent = true;
                     break;
                  }
               }
               faceSet.clear();
               localBlock->set_extent(iLocalVert, jLocalVert, kLocalVert,
                                      atDomainExtent);
               for ( int iNbrBlk = 0; iNbrBlk != numNbrVert; ++iNbrBlk )
                  nbrList[iNbrBlk]->set_extent(nbrIndices[iNbrBlk].i,
                                               nbrIndices[iNbrBlk].j,
                                               nbrIndices[iNbrBlk].k,
                                               atDomainExtent);
            }

//--Only compute transformations for blocks unique to the vertex

            Sort::intro(numNbrVert, nbrList);
            numNbrVert = erase_common(nbrList, nbrList + numNbrVert,
                                      nbrEdgeList, nbrEdgeList + numNbrAllEdge);

            for ( int iNbrBlk = 0; iNbrBlk != numNbrVert; ++iNbrBlk ) {
               Block *donorBlock = nbrList[iNbrBlk];

               // Check if this vertex neighbour has yet been set
               if ( !localBlock->have_neighbour(
                       iLocalVert, jLocalVert, kLocalVert,
                       donorBlock->label) ) {

                  // Coordinates of the donor vertex described by ijk inidices
                  const int mDonorVert = donorBlock->find_vertex_index(vert);
                  int iDonorVert, jDonorVert, kDonorVert;
                  Block::mvtoijk(mDonorVert,
                                 iDonorVert, jDonorVert, kDonorVert);

//--Describe the transformation

                  // Compact transformation matrices describing the orientation
                  // of the local block with respect to the reference frame.
                  // This describes the orientation of the donor block if the
                  // vertex index is taken as [7-mv]
                  static const int CTMRefToLocalVert[8][3] = {
                     { 1,  2,  3}, {-1,  3,  2}, { 1,  3, -2}, {-1, -2,  3},
                     { 1, -3,  2}, {-1,  2, -3}, { 1, -2, -3}, {-1, -3, -2}
                  };

                  // Adjustments to the donor block orientation based on
                  // rotation around the shared vertex
                  static const int rotCTMVert[3][3] = {
                     {1, 2, 3}, {3, 1, 2}, {2, 3, 1}
                  };

                  // Vectors in the reference i, j, and k directions for the
                  // local block.  This is the reference for aligning the
                  // vectors.
                  double vecSetLocal[3][3];
                  localBlock->get_vector(
                     mLocalVert, std::abs(CTMRefToLocalVert[mLocalVert][0]) - 1,
                     vecSetLocal[0]);
                  localBlock->get_vector(
                     mLocalVert, std::abs(CTMRefToLocalVert[mLocalVert][1]) - 1,
                     vecSetLocal[1]);
                  localBlock->get_vector(
                     mLocalVert, std::abs(CTMRefToLocalVert[mLocalVert][2]) - 1,
                     vecSetLocal[2]);

                  // Vectors in the reference i, j, and k directions for the
                  // donor block
                  double vecSetDonor[3][3];
                  donorBlock->get_vector(
                     mDonorVert,
                     std::abs(CTMRefToLocalVert[7-mDonorVert][0]) - 1,
                     vecSetDonor[0]);
                  donorBlock->get_vector(
                     mDonorVert,
                     std::abs(CTMRefToLocalVert[7-mDonorVert][1]) - 1,
                     vecSetDonor[1]);
                  donorBlock->get_vector(
                     mDonorVert,
                     std::abs(CTMRefToLocalVert[7-mDonorVert][2]) - 1,
                     vecSetDonor[2]);

                  // Weights of a dot-product between different orientations of
                  // the vectors.
                  double wtRot[3];
                  wtRot[0] = set_dot(vecSetLocal, vecSetDonor, 0);
                  wtRot[1] = set_dot(vecSetLocal, vecSetDonor, 1);
                  wtRot[2] = set_dot(vecSetLocal, vecSetDonor, 2);

                  // The minimum weight gives the closest alignment.  Specify
                  // the amount of rotation around the shared vertex
                  int iRot = 0;
                  if ( wtRot[1] < wtRot[0] ) iRot = 1;
                  if ( wtRot[2] < wtRot[iRot] ) iRot = 2;
                  int rotIndex;         // Index into CTMRefToLocalVert (but for
                                        // donor blocks) as modified by rotation
                  int rotSign;          // Sign change by rotation.  The sign
                                        // is not changed by vertex rotation and
                                        // is therefore not used.
                  int ctm[6];           // Compact transformation matrix

                  // Transformations between the local and reference
                  ctm[0] = CTMRefToLocalVert[mLocalVert][0];
                  ctm[1] = CTMRefToLocalVert[mLocalVert][1];
                  ctm[2] = CTMRefToLocalVert[mLocalVert][2];
                  ctm[3] = 0;
                  ctm[4] = 0;
                  ctm[5] = 0;

                  CoordTransform TMRefToLocal(ctm);
                  CoordTransform TMLocalToRef(TMRefToLocal);
                  TMLocalToRef.transpose();

                  // Transformations between the donor and reference
                  rotate_value(rotCTMVert[iRot][0], rotSign, rotIndex);
                  ctm[0] = CTMRefToLocalVert[7-mDonorVert][rotIndex];
                  rotate_value(rotCTMVert[iRot][1], rotSign, rotIndex);
                  ctm[1] = CTMRefToLocalVert[7-mDonorVert][rotIndex];
                  rotate_value(rotCTMVert[iRot][2], rotSign, rotIndex);
                  ctm[2] = CTMRefToLocalVert[7-mDonorVert][rotIndex];

                  CoordTransform TMRefToDonor(ctm);
                  CoordTransform TMLocalToDonor(TMRefToDonor, TMLocalToRef);

//--Describe the offset.

                  // Offset from local block vertex 0 to reference frame 0.
                  static const int OffsetLocalToRefVert[8][3] = {
                     { 0,  0,  0}, {-1,  0,  0}, {0,  0, -1}, {-1, -1,  0},
                     { 0, -1,  0}, {-1,  0, -1}, {0, -1, -1}, {-1, -1, -1}
                  };

                  // Offset from donor block vertex 0 to reference frame 0.
                  static const int OffsetDonorToRefVert[8][3] = {
                     { 0,  0,  0}, { 1,  0,  0}, {0,  1,  0}, { 1,  0,  1},
                     { 0,  0,  1}, { 1,  1,  0}, {0,  1,  1}, { 1,  1,  1}
                  };

                  // rotOffsetVert is the same as rotCTMVert

                  // Offset for local block.
                  int localOffset[3];
                  localOffset[0] = OffsetLocalToRefVert[mLocalVert][0];
                  localOffset[1] = OffsetLocalToRefVert[mLocalVert][1];
                  localOffset[2] = OffsetLocalToRefVert[mLocalVert][2];

                  // The localOffset, defined in the reference frame, needs to
                  // be transformed into the local frame, just so it can be
                  // multiplied by the block dimensions.  Then it needs to be
                  // tranformed to the donor frame.
                  TMRefToLocal.transform(localOffset);  // To local frame
                  localOffset[0] *= localBlock->dim[0];
                  localOffset[1] *= localBlock->dim[1];
                  localOffset[2] *= localBlock->dim[2];
                  // Considerations for ghost cells
                  localOffset[0] += numGhostLayer;
                  localOffset[1] += numGhostLayer;
                  localOffset[2] += numGhostLayer;
                  TMLocalToDonor.transform(localOffset);  // Put in donor frame

                  // Offset for donor block.  The donor offset is rotated by the
                  // pattern rotCTMVert
                  int donorOffset[3];
                  rotate_value(rotCTMVert[iRot][0], rotSign, rotIndex);
                  donorOffset[0] = OffsetDonorToRefVert[mDonorVert][rotIndex];
                  rotate_value(rotCTMVert[iRot][1], rotSign, rotIndex);
                  donorOffset[1] = OffsetDonorToRefVert[mDonorVert][rotIndex];
                  rotate_value(rotCTMVert[iRot][2], rotSign, rotIndex);
                  donorOffset[2] = OffsetDonorToRefVert[mDonorVert][rotIndex];

                  // The donorOffset, defined in the reference frame, needs to
                  // be transformed into the donor frame and multiplied by the
                  // block dimensions.
                  TMRefToDonor.transform(donorOffset);  // To donor frame
                  donorOffset[0] *= donorBlock->dim[0];
                  donorOffset[1] *= donorBlock->dim[1];
                  donorOffset[2] *= donorBlock->dim[2];
                  // Considerations for ghost cells
                  donorOffset[0] += numGhostLayer;
                  donorOffset[1] += numGhostLayer;
                  donorOffset[2] += numGhostLayer;

                  // For cell-center transformations, add a jump across the
                  // boundary
                  if ( cellType == CellCenter ) {
                     donorOffset[0] += iDonorVert;
                     donorOffset[1] += jDonorVert;
                     donorOffset[2] += kDonorVert;
                  }

                  // Subtract the local_to_ref offset from the donor_to_ref
                  // offset to get the total donor_to_local offset
                  int totalOffset[3];
                  totalOffset[0] = donorOffset[0] - localOffset[0];
                  totalOffset[1] = donorOffset[1] - localOffset[1];
                  totalOffset[2] = donorOffset[2] - localOffset[2];
                  TMLocalToDonor.set_offset(totalOffset);

                  // Get the reverse transformation
                  CoordTransform TMDonorToLocal(TMLocalToDonor);
                  TMDonorToLocal.reverse();

//--Set the neighbour information

                  // Set the donor information as seen from the local block
                  NeighbourData *nbrData = nbrDataMem->get();
                  nbrData->label = donorBlock->label;
                  nbrData->boundary[0] = iDonorVert;
                  nbrData->boundary[1] = jDonorVert;
                  nbrData->boundary[2] = kDonorVert;
                  TMLocalToDonor.compact(nbrData->ctm);
                  localBlock->set_neighbour_data(iLocalVert, jLocalVert,
                                                 kLocalVert, nbrData);

                  // Set the local information as seen from the donor block
                  nbrData = nbrDataMem->get();
                  nbrData->label = localBlock->label;
                  nbrData->boundary[0] = iLocalVert;
                  nbrData->boundary[1] = jLocalVert;
                  nbrData->boundary[2] = kLocalVert;
                  TMDonorToLocal.compact(nbrData->ctm);
                  donorBlock->set_neighbour_data(iDonorVert, jDonorVert,
                                                 kDonorVert, nbrData);
               
               }  // End if this neighbour data not yet set
            }  // Loop over neighbour blocks across this vertex
         }  // End of loop over local vertices
         
      }  // End if block exists
   }  // End of loop over blocks

   delete[] nbrIndices;
   delete[] nbrList;
   delete[] nbrFaceList;
   delete[] nbrEdgeList;

}


/*==============================================================================
 *
 * Routine num_neighbour_block
 *
 * Purpose
 * =======
 *
 *   Returns the number of neighbour blocks across a boundary element.
 *
 * I/O
 * ===
 *
 *   blockIndex         - (I) Index of local block
 *   iloc, jloc, kloc   - (I) Index of boundary element across which to find
 *                            number of neighbours
 *   return             - number of neighbour blocks
 *
 * Notes
 * =====
 *
 *   - if zero is returned, consider calling domain_extent to see if there are
 *     zero neighbours because this boundary element is at a domain extent.
 *     Otherwise, it is because of the unstructured connectivity between the
 *     blocks.
 *
 *============================================================================*/

int BlkC::BlockConnectivity::num_neighbour_block(const int blockIndex,
                                                 const int iloc,
                                                 const int jloc,
                                                 const int kloc)
{

   if ( !computeCalled ) compute_neighbours();
   Block *const blkList = static_cast<Block*>(_blkList);
   return blkList[blockIndex].num_neighbour(iloc, jloc, kloc);

}


/*==============================================================================
 *
 * Routine neighbour_block
 *
 * Purpose
 * =======
 *
 *   Returns data for a neighbour block across a boundary element
 *
 * I/O
 * ===
 *
 *   blockIndex         - (I) Index of local block
 *   iloc, jloc, kloc   - (I) Index of element across which to find neighbour
 *   nbrloc             - (I) Neighbour to find.  Start counting from 0.
 *   blockIndexNbr      - (O) Index of the neighbour block
 *   ilocNbr, jlocNbr, klocNbr
 *                      - (O) Index of the shared boundary element in the
 *                            neighbour block.
 *   ctm                - (O) A rank 1 array of of dimension 6
 *                            [0] to [2] the compact transformation
 *                            [3] to [5] the offset
 *   return             - 0 - success
 *                        1 - requested neighbour not found
 *
 *============================================================================*/

int BlkC::BlockConnectivity::neighbour_block(const int blockIndex,
                                             const int iloc,
                                             const int jloc,
                                             const int kloc,
                                             const int nbrloc,
                                             int &blockIndexNbr,
                                             int &ilocNbr,
                                             int &jlocNbr,
                                             int &klocNbr,
                                             int *const ctm)
{

   if ( !computeCalled ) compute_neighbours();
   Block *const blkList = static_cast<Block*>(_blkList);
   return blkList[blockIndex].neighbour_data(iloc, jloc, kloc,
                                             nbrloc, blockIndexNbr,
                                             ilocNbr, jlocNbr, klocNbr,
                                             ctm);

}


/*==============================================================================
 *
 * Routine at_domain_extent
 *
 * Purpose
 * =======
 *
 *   Checks if the boundary element of a block is at the domain extent.
 *
 * I/O
 * ===
 *
 *   blockIndex         - (I) Index of local block
 *   iloc, jloc, kloc   - (I) Index of boundary element across which to find
 *                            number of neighbours
 *   return             - T - at domain extent
 *                        F - internal
 *
 * Notes
 * =====
 *
 *   - this is useful for determining if block has no neighbours across a
 *     given boundary element (vertex, edge, or face) because of unstructured
 *     connectivity or because it is at the edge of the domain.
 *
 *============================================================================*/

bool BlkC::BlockConnectivity::at_domain_extent(const int blockIndex,
                                               const int iloc,
                                               const int jloc,
                                               const int kloc)
{

   if ( !computeCalled ) compute_neighbours();
   Block *const blkList = static_cast<Block*>(_blkList);
   return blkList[blockIndex].at_extent(iloc, jloc, kloc);

}


/*******************************************************************************
 *
 * Member routines for debugging class BlockConnectivity
 *
 ******************************************************************************/

#ifdef BLKC_ENABLE_DEBUGGING


/*==============================================================================
 *
 * Routine dbg_save_input
 *
 * Purpose
 * =======
 *
 *   Saves information about the blocks that have been loaded into the database
 *
 * I/O
 * ===
 *
 *   fileName           - (I) Name of file for saving block information
 *
 * Notes
 * =====
 *
 *   - Call this routine before any blocks have been added
 *   - The destructor will close the file
 *
 *============================================================================*/

void BlkC::BlockConnectivity::dbg_save_input(const char *const fileName)
{
   fout.open(fileName);
   fout << sizeBlkList << ' ' << cellType << ' ' << numGhostLayer << std::endl;
}

#endif  // BLKC_ENABLE_DEBUGGING
