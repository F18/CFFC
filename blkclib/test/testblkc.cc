
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

#include <iostream>
#include <iomanip>
#include "block_connectivity.h"
#include "coord_transform.h"


/*==============================================================================
 * A simple vector class
 *============================================================================*/

template <typename T>
class Vec
{
  private:
   T _vec[3];

  public:
   Vec() 
   {
      _vec[0] = 0;
      _vec[1] = 0;
      _vec[2] = 0;
   }
   Vec(const T *const vec)
   {
      vec[0] = vec[0];
      vec[1] = vec[1];
      vec[2] = vec[2];
   }
   Vec(const T v0, const T v1, const T v2)
   {
      _vec[0] = v0;
      _vec[1] = v1;
      _vec[2] = v2;
   }
   // Use sythesized copy, assignment, and destructor

   T x() const { return _vec[0]; }
   T &x() { return _vec[0]; }

   T y() const { return _vec[1]; }
   T &y() { return _vec[1]; }

   T z() const { return _vec[2]; }
   T &z() { return _vec[2]; }

   const T *get_vector() const { return _vec; }
   T *set_vector() { return _vec; }
};

template <typename T>
inline Vec<T> operator-(const Vec<T> &a, const Vec<T> &b)
{
   return Vec<T>(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
}


/*==============================================================================
 * A simple representation of a grid block
 *============================================================================*/

class Block
{
  public:
   Vec<double> coord[2][2][2];
   // Use sythesized constructor, copy, assignment, and destructor
};


/*******************************************************************************
 *
 * Routine main
 *
 * Purpose
 * =======
 *
 *   Tests the BlockConnectivity class across all vertices, edges, and faces,
 *   and every orientation.  Basically, a center block, shown shaded in the
 *   figure below is wrapped by 26 neighbour blocks.  The neighbour blocks
 *   all have an orientation aligned with the reference orientation (shown by
 *   (i, j, k) in the figure).  The center block is rotated through all 24
 *   possible orientations.  The data about each neighbour block as returned by
 *   the BlockConnectivity class for the center block is examined to see if it
 *   is correct.  The reverse information (data about the center as returned by
 *   the BlockConnectivity class for each neighbour) is also checked.  Results
 *   are reported to standard output.
 *
 *                     ^ j
 *                     |
 *                     |
 *                     |
 *                     o----------o----------o----------o
 *                    /          /          /          /|
 *                   /          /          /          / |
 *                  /          /          /          /  |
 *                 o----------o----------o----------o   |
 *                /          /|          |          |   |
 *               /          / |          |          |   o
 *              /          /  |          |          |  /|
 *             o----------o   |          |          | / |
 *            /          /|   |          |          |/  |
 *           /          / |   o----------o----------o   |
 *          /          /  |  ////////////|          |   |
 *         o----------o   | ////////////||          |   o
 *         |          |   |////////////|||          |  /|
 *         |          |   o----------o||||          | / |
 *         |          |  /|----------|||||          |/  |
 *         |          | / |----------||||o----------o   |
 *         |          |/  |----------|||/          /|   |       i
 *         o----------o   |----------||/          / |   o ------>
 *         |          |   |----------|/          /  |  /
 *         |          |   o----------o----------o   | /
 *         |          |  /          /          /|   |/
 *         |          | /          /          / |   o
 *         |          |/          /          /  |  /
 *         o----------o----------o----------o   | /
 *         |          |          |          |   |/
 *         |          |          |          |   o
 *         |          |          |          |  /
 *         |          |          |          | /
 *         |          |          |          |/
 *         o----------o----------o----------o
 *        /
 *       /
 *      / k
 *     `-
 *
 ******************************************************************************/

int main()
{

   // Possible orientations of the center block
   static const int orientation[24][3] = {
      {-1,  3,  2}, {-1,  2, -3}, {-1, -3, -2}, {-1, -2,  3},
      { 1,  2,  3}, { 1,  3, -2}, { 1, -2, -3}, { 1, -3,  2},
      {-2,  1,  3}, {-2,  3, -1}, {-2, -1, -3}, {-2, -3,  1},
      { 2,  3,  1}, { 2,  1, -3}, { 2, -3, -1}, { 2, -1,  3},
      {-3,  2,  1}, {-3,  1, -2}, {-3, -2, -1}, {-3, -1,  2},
      { 3,  1,  2}, { 3,  2, -1}, { 3, -1, -2}, { 3, -2,  1}
   };
   Block block[3][3][3];

//--Define the outer blocks

   for ( int iBlk = 0; iBlk != 3; ++iBlk ) {
      for ( int jBlk = 0; jBlk != 3; ++jBlk ) {
         for ( int kBlk = 0; kBlk != 3; ++kBlk ) {
            if ( iBlk*jBlk*kBlk != 1 ) {

               for ( int i = 0; i != 2; ++i ) {
                  for ( int j = 0; j != 2; ++j ) {
                     for ( int k = 0; k != 2; ++k ) {
                        block[iBlk][jBlk][kBlk].coord[i][j][k].x() =
                           -1.5 + iBlk + i;
                        block[iBlk][jBlk][kBlk].coord[i][j][k].y() =
                           -1.5 + jBlk + j;
                        block[iBlk][jBlk][kBlk].coord[i][j][k].z() =
                           -1.5 + kBlk + k;
                     }
                  }
               }

            }
         }
      }
   }


/*==============================================================================
 * Loop over all orientations
 *============================================================================*/

   bool globalError = false;
   for ( int mOr = 0; mOr != 24; ++mOr ) {
      const int iOr = orientation[mOr][0];
      const int jOr = orientation[mOr][1];
      const int kOr = orientation[mOr][2];

//--Set the coordinates of the local block

      int ctm[6];
      ctm[0] = iOr;
      ctm[1] = jOr;
      ctm[2] = kOr;
      ctm[std::abs(iOr)+2] = ( iOr < 0 ) ? 1 : 0;
      ctm[std::abs(jOr)+2] = ( jOr < 0 ) ? 1 : 0;
      ctm[std::abs(kOr)+2] = ( kOr < 0 ) ? 1 : 0;
      CoordTransform TMRefToCen(ctm);

      for ( int i = 0; i != 2; ++i ) {
         for ( int j = 0; j != 2; ++j ) {
            for ( int k = 0; k != 2; ++k ) {
               int index[3];
               index[0] = i;
               index[1] = j;
               index[2] = k;
               TMRefToCen.transform(index);
               block[1][1][1].coord[index[0]][index[1]][index[2]].x() =
                  -0.5 + i;
               block[1][1][1].coord[index[0]][index[1]][index[2]].y() =
                  -0.5 + j;
               block[1][1][1].coord[index[0]][index[1]][index[2]].z() =
                  -0.5 + k;
            }
         }
      }

/*--------------------------------------------------------------------*
 * Test the block connectivity class
 *--------------------------------------------------------------------*/

      BlkC::BlockConnectivity blockConn(27, BlkC::CellCenter, 2, 6);
      BlkC::VertexPack vp;
      int mBlock = 0;

#ifdef BLKC_ENABLE_DEBUGGING
      if ( mOr == 0 ) blockConn.dbg_save_input("test_input.dat");
#endif

//--Add the center block
      
      {
         Block &wBlock = block[1][1][1];
         for ( int i = 0; i != 2; ++i ) {
            for ( int j = 0; j != 2; ++j ) {
               for ( int k = 0; k != 2; ++k ) {
                  BlkC::ILoc_t iloc = ( i ) ? BlkC::IMax : BlkC::IMin;
                  BlkC::JLoc_t jloc = ( j ) ? BlkC::JMax : BlkC::JMin;
                  BlkC::KLoc_t kloc = ( k ) ? BlkC::KMax : BlkC::KMin;
                  // Set coordinate
                  vp.set_coord(iloc, jloc, kloc,
                               wBlock.coord[i][j][k].get_vector());
                  // Set vectors along each edge
                  const int iOff = ( i ) ? 0 : 1;
                  const int jOff = ( j ) ? 0 : 1;
                  const int kOff = ( k ) ? 0 : 1;
                  Vec<double> iVec, jVec, kVec;
                  iVec = wBlock.coord[iOff][j][k] - wBlock.coord[i][j][k];
                  jVec = wBlock.coord[i][jOff][k] - wBlock.coord[i][j][k];
                  kVec = wBlock.coord[i][j][kOff] - wBlock.coord[i][j][k];
                  vp.set_vector(iloc, jloc, kloc, 'i', iVec.get_vector());
                  vp.set_vector(iloc, jloc, kloc, 'j', jVec.get_vector());
                  vp.set_vector(iloc, jloc, kloc, 'k', kVec.get_vector());
               }
            }
         }
      }

      blockConn.add_block(mBlock++, 5, 5, 5, vp);

//--Add all the other blocks

      for ( int iBlk = 0; iBlk != 3; ++iBlk ) {
         for ( int jBlk = 0; jBlk != 3; ++jBlk ) {
            for ( int kBlk = 0; kBlk != 3; ++kBlk ) {
               if ( iBlk*jBlk*kBlk != 1 ) {

                  Block &wBlock = block[iBlk][jBlk][kBlk];
                  for ( int i = 0; i != 2; ++i ) {
                     for ( int j = 0; j != 2; ++j ) {
                        for ( int k = 0; k != 2; ++k ) {
                           BlkC::ILoc_t iloc = ( i ) ? BlkC::IMax : BlkC::IMin;
                           BlkC::JLoc_t jloc = ( j ) ? BlkC::JMax : BlkC::JMin;
                           BlkC::KLoc_t kloc = ( k ) ? BlkC::KMax : BlkC::KMin;
                           // Set coordinate
                           vp.set_coord(iloc, jloc, kloc,
                                        wBlock.coord[i][j][k].get_vector());
                           // Set vectors along each edge
                           const int iOff = ( i ) ? 0 : 1;
                           const int jOff = ( j ) ? 0 : 1;
                           const int kOff = ( k ) ? 0 : 1;
                           Vec<double> iVec, jVec, kVec;
                           iVec = wBlock.coord[iOff][j][k] -
                              wBlock.coord[i][j][k];
                           jVec = wBlock.coord[i][jOff][k] -
                              wBlock.coord[i][j][k];
                           kVec = wBlock.coord[i][j][kOff] -
                              wBlock.coord[i][j][k];
                           vp.set_vector(iloc, jloc, kloc, 'i',
                                         iVec.get_vector());
                           vp.set_vector(iloc, jloc, kloc, 'j',
                                         jVec.get_vector());
                           vp.set_vector(iloc, jloc, kloc, 'k',
                                         kVec.get_vector());
                        }
                     }
                  }

                  blockConn.add_block(mBlock++, 5, 5, 5, vp);

               }
            }
         }
      }

//--Check extents

      if ( mOr == 0 ) {

         std::cout << "Checking extents ...";
         bool errorExtents = false;
         const char* LBLextent[] = { "Interior", "Domain extent" };

         // Check the center block (should have no boundary elements at the
         // extents of the domain.
         for ( int i = -1; i != 2; ++i ) {
            for ( int j = -1; j != 2; ++j ) {
               for ( int k = -1; k != 2; ++k ) {
                  if ( std::abs(i) + std::abs(j) + std::abs(k) ) {
                     if ( blockConn.at_domain_extent(0, i, j, k) ) {
                        errorExtents = true;
                        std::cout << "\n  Error: extent not correct for center "
                           "block at boundary ("
                                  << std::setw(2) << i << ", "
                                  << std::setw(2) << j << ", "
                                  << std::setw(2) << k << ')';
                        std::cout << "\n    Expected: " << LBLextent[0];
                        std::cout << "\n    Found   : " << LBLextent[1];
                     }
                  }
               }
            }
         }

         // Check all other blocks
         mBlock = 1;
         for ( int iBlk = 0; iBlk != 3; ++iBlk ) {
            for ( int jBlk = 0; jBlk != 3; ++jBlk ) {
               for ( int kBlk = 0; kBlk != 3; ++kBlk ) {
                  if ( iBlk*jBlk*kBlk != 1 ) {

                     // Indicate expected extents of the block (99 means the
                     // block is interior on all boundaries in this compuational
                     // index)
                     const int iExt = ( iBlk == 1 ) ? 99 : iBlk - 1;
                     const int jExt = ( jBlk == 1 ) ? 99 : jBlk - 1;
                     const int kExt = ( kBlk == 1 ) ? 99 : kBlk - 1;
                     for ( int i = -1; i != 2; ++i ) {
                        for ( int j = -1; j != 2; ++j ) {
                           for ( int k = -1; k != 2; ++k ) {
                              if ( std::abs(i) + std::abs(j) + std::abs(k) ) {
                                 bool extentFound =
                                    blockConn.at_domain_extent(mBlock, i, j, k);
                                 bool extentExpected = false;
                                 if ( i == iExt || j == jExt || k == kExt )
                                    extentExpected = true;
                                 if ( extentFound != extentExpected ) {
                                    errorExtents = true;
                                    std::cout << "\n  Error: extent not "
                                       "correct for block ("
                                              << std::setw(1) << iBlk << ", "
                                              << std::setw(1) << jBlk << ", "
                                              << std::setw(1) << kBlk << ") at "
                                       "boundary ("
                                              << std::setw(2) << i << ", "
                                              << std::setw(2) << j << ", "
                                              << std::setw(2) << k << ')';
                                    std::cout << "\n    Expected: "
                                              << LBLextent[extentExpected];
                                    std::cout << "\n    Found   : "
                                              << LBLextent[extentFound];
                                 }
                              }
                           }
                        }
                     }
                     ++mBlock;

                  }
               }
            }
         }

         if ( errorExtents ) {
            globalError = true;
            std::cout << "\n                                                   "
               "                   FAILED\n";
         }
         else {
            std::cout << "                                                  PAS"
               "SED\n";
         }

      }

//--Check orientations

      std::cout << "Checking orientation ("
                << std::setw(2) << iOr << ", "
                << std::setw(2) << jOr << ", "
                << std::setw(2) << kOr << ") ...";
      bool errorThisOrientation = false;
      mBlock = 1;
      TMRefToCen.clear_offset();
      for ( int iBlk = 0; iBlk != 3; ++iBlk ) {
         for ( int jBlk = 0; jBlk != 3; ++jBlk ) {
            for ( int kBlk = 0; kBlk != 3; ++kBlk ) {
               if ( iBlk*jBlk*kBlk != 1 ) {

                  bool errorThisNbrBlk = false;
                  // Get location of center block with respect to this block
                  int cenBlkLoc[3];
                  cenBlkLoc[0] = 1 - iBlk;
                  cenBlkLoc[1] = 1 - jBlk;
                  cenBlkLoc[2] = 1 - kBlk;
                  // Get this block location with respect to the center block
                  int blkLocWRTCen[3];
                  blkLocWRTCen[0] = iBlk - 1;
                  blkLocWRTCen[1] = jBlk - 1;
                  blkLocWRTCen[2] = kBlk - 1;
                  TMRefToCen.transform(blkLocWRTCen);
                  // Query information about this neighbour from the center
                  int indexNbr;
                  int iLocNbr, jLocNbr, kLocNbr;
                  int cctm[6];
                  blockConn.neighbour_block(0, blkLocWRTCen[0], blkLocWRTCen[1],
                                            blkLocWRTCen[2], 0, indexNbr,
                                            iLocNbr, jLocNbr, kLocNbr, cctm);
                  if ( indexNbr != mBlock ) {
                     errorThisNbrBlk = true;
                     std::cout << "\n  Error: neighbour index is not correct";
                     std::cout << "\n    Expected: " << mBlock << std::endl;
                     std::cout << "\n    Found   : " << indexNbr << std::endl;
                  }
                  if ( iLocNbr != cenBlkLoc[0] ||
                       jLocNbr != cenBlkLoc[1] ||
                       kLocNbr != cenBlkLoc[2] ) {
                     errorThisNbrBlk = true;
                     std::cout << "\n  Error: neighbour boundary is not "
                        "correct";
                     std::cout << "\n    Expected: ("
                               << std::setw(2) << cenBlkLoc[0] << ", "
                               << std::setw(2) << cenBlkLoc[1] << ", "
                               << std::setw(2) << cenBlkLoc[2] << ')';
                     std::cout << "\n    Found   : ("
                               << std::setw(2) << iLocNbr << ", "
                               << std::setw(2) << jLocNbr << ", "
                               << std::setw(2) << kLocNbr << ')';
                  }

                  // We want to find index (3, 4, 2) in each neighbour block.
                  // This index in the center block is given for the '3' sets
                  // of neighbour blocks in each of the i, j, and k reference
                  // directions.  The index in the center block depends on which
                  // face aligns with the imax, jmax, and kmax reference faces.
                  static const int iSetMax[3] = {-2,  3,  8};
                  static const int iSetMin[3] = {10,  5,  0};
                  static const int jSetMax[3] = {-1,  4,  9};
                  static const int jSetMin[3] = { 9,  4, -1};
                  static const int kSetMax[3] = {-3,  2,  7};
                  static const int kSetMin[3] = {11,  6,  1};

                  int nIndex[3];        // Index in the neighbour block
                  int cIndex[3];        // Index in the center block
                  const int iSetIndex = std::abs(iOr) - 1;
                  const int jSetIndex = std::abs(jOr) - 1;
                  const int kSetIndex = std::abs(kOr) - 1;
                  cIndex[iSetIndex] =
                     ( iOr < 0 ) ? iSetMin[iBlk] : iSetMax[iBlk];
                  cIndex[jSetIndex] =
                     ( jOr < 0 ) ? jSetMin[jBlk] : jSetMax[jBlk];
                  cIndex[kSetIndex] =
                     ( kOr < 0 ) ? kSetMin[kBlk] : kSetMax[kBlk];
                  nIndex[0] = cIndex[0];
                  nIndex[1] = cIndex[1];
                  nIndex[2] = cIndex[2];
                  CoordTransform TMCenToNbr(cctm);
                  TMCenToNbr.transform(nIndex);
                  if ( nIndex[0] != 3 || nIndex[1] != 4 || nIndex[2] != 2 ) {
                     errorThisNbrBlk = true;
                     std::cout << "\n  Error: Index (3, 4, 2) not computed by "
                        "transformation";
                     std::cout << "\n    Instead found: ("
                               << std::setw(2) << nIndex[0] << ", "
                               << std::setw(2) << nIndex[1] << ", "
                               << std::setw(2) << nIndex[2] << ')';
                  }

                  // Check the reverse transformation
                  int nctm[3];
                  blockConn.neighbour_block(mBlock, cenBlkLoc[0], cenBlkLoc[1],
                                            cenBlkLoc[2], 0, indexNbr, iLocNbr,
                                            jLocNbr, kLocNbr, nctm);
                  if ( nctm[0] != iOr || nctm[1] != jOr || nctm[2] != kOr ) {
                     errorThisNbrBlk = true;
                     std::cout << "\n  Error: Neighbour transformation matrix "
                        "is not the same as reference";
                     std::cout << "\n    Expected: ("
                               << std::setw(2) << iOr << ", "
                               << std::setw(2) << jOr << ", "
                               << std::setw(2) << kOr << ')';
                     std::cout << "\n    Found   : ("
                               << std::setw(2) << nctm[0] << ", "
                               << std::setw(2) << nctm[1] << ", "
                               << std::setw(2) << nctm[2] << ')';
                  }

                  CoordTransform TMNbrToCen(nctm);
                  TMNbrToCen.transform(nIndex);
                  if ( nIndex[0] != cIndex[0] || nIndex[1] != cIndex[1] ||
                       nIndex[2] != cIndex[2] ) {
                     errorThisNbrBlk = true;
                     std::cout << "\n  Error: Reverse transformation incorrect";
                     std::cout << "\n    Expected: ("
                               << std::setw(2) << cIndex[0] << ", "
                               << std::setw(2) << cIndex[1] << ", "
                               << std::setw(2) << cIndex[2] << ')';
                     std::cout << "\n    Found   : ("
                               << std::setw(2) << nIndex[0] << ", "
                               << std::setw(2) << nIndex[1] << ", "
                               << std::setw(2) << nIndex[2] << ')';
                  }

                  if ( errorThisNbrBlk ) {
                     errorThisOrientation = true;
                     std::cout << "\n    ... at neighbour block ("
                               << std::setw(2) << iBlk - 1 << ", "
                               << std::setw(2) << jBlk - 1 << ", "
                               << std::setw(2) << kBlk - 1 << ')';
                  }
                  ++mBlock;

               }
            }
         }
      }

      if ( errorThisOrientation ) {
         globalError = true;
         std::cout << "\n                                                      "
            "                FAILED\n";
      }
      else {
         std::cout << "                                 PASSED\n";
      }

   }  // End loop over orientations

   if ( !globalError ) std::cout << "All checks passed!\n";

}
