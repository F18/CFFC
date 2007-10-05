#ifndef _COORD_TRANSFORM_INCLUDED
#define _COORD_TRANSFORM_INCLUDED

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

//--Forward declarations

namespace BlkC
{
   class BlockConnectivity;
}


/*******************************************************************************
 *
 * class CoordTransform
 *
 * Purpose
 * =======
 *
 *   Provides a matrix and functions for transforming indices from the local
 *   coordinate system of one block to another.
 *
 * Constructors
 * ============
 *
 *   CoordTransform()  -- default constructor sets matrix to identity and offset
 *                        to zero
 *   CoordTransform(const int *const ctm)
 *                     -- constructs matrix and offset from compact storage
 *     ctm                (I) compact transformation matrix and offset.  Rank 1
 *                            array with dimension (6).
 *
 * Member Functions
 * ================
 *
 *   void expand(const int *const ctm)
 *                     -- expands a compact transformation matrix into this
 *                        matrix and offset
 *     ctm                (I) compact transformation matrix and offset.  Rank 1
 *                            array with dimension (6).
 *
 *   void compact(int *const ctm) const
 *                     -- stores this matrix and offset into compact storage
 *     ctm                (O) compact transformation matrix and offset.  Rank 1
 *                            array with dimension (6).
 *
 *   void transform(int *const index) const
 *                     -- transforms index
 *     index              (I) rank 1 array with dimension (3) describing index
 *                            in this coordinate system.
 *                        (O) index in neighbour coordinate system
 *
 *   void transform(int &i, int &j, int &k) const
 *                     -- transforms index
 *     i, j, k            (I) indices in this coordinate system.
 *                        (O) indices in neighbour coordinate system.
 *
 *   void transpose()  -- provides transpose of this matrix.  The transpose
 *                        allows for transformations in the opposite direction
 *
 *   void reverse()    -- provides the reverse coordinate transformation.
 *                        Whereas transform() simply modifies the matrix,
 *                        reverse() also modifies the offset.
 *
 *   void clear_offset()
 *                     -- sets the offset to 0.  This is useful for transforming
 *                        coordinate frames instead of indices.
 *
 * Member Operators
 * ================
 *
 *   int &operator()(const int i, const int j) 
 *                     -- indexes the matrix
 *     i, j               (I) indices
 *     return             (O) reference to matrix element
 *
 *   const int &operator()(const int i, const int j) const
 *                     -- indexes the matrix
 *     i, j               (I) indices
 *     return             (O) reference to matrix element
 *
 *   void operator=(const int val)
 *                     -- assigns value to entire matrix
 *     val                (I) value to assign
 *
 * Notes
 * =====
 *
 *   - See the CGNS docs:
 *     <http://www.grc.nasa.gov/WWW/cgns/sids/cnct.html#GridConnectivity1to1>
 *     for details on the transformation matrix.
 *   - In CGNS, the elements of the compact storage matrix has a range of 1 to
 *     3.  Here, the compact storage matrix has a range of 0 to 2.
 *   - Row storage is used for efficient matrix-vector multiplication
 *     +-     -+
 *     | 0 1 2 |
 *     | 3 4 5 |
 *     | 6 7 8 |
 *     +-     -+
 *
 ******************************************************************************/

class CoordTransform
{


/*==============================================================================
 * Data
 *============================================================================*/

  private:
   int _tm[9];                          // Transformation matrix
   int _os[3];                          // Offset.  Same as begin2 in CGNS docs
                                        // if begin1 = (0, 0, 0)


/*==============================================================================
 * Public member functions
 *============================================================================*/

  public:

//--Constructors

   // Default
   CoordTransform()
   {
      _tm[0] = 1;
      _tm[1] = 0;
      _tm[2] = 0;
      _tm[3] = 0;
      _tm[4] = 1;
      _tm[5] = 0;
      _tm[6] = 0;
      _tm[7] = 0;
      _tm[8] = 1;
      clear_offset();
   }

   // From compact storage
   CoordTransform(const int *const ctm)
   {
      expand(ctm);
   }

   // Use synthesized copy, assignment

//--Use synthesized destructor

//--Expand from compact storage

   void expand(const int *const ctm)
   {
      this->operator=(0);
      this->operator()(std::abs(ctm[0]) - 1, 0) = sign(ctm[0]);
      this->operator()(std::abs(ctm[1]) - 1, 1) = sign(ctm[1]);
      this->operator()(std::abs(ctm[2]) - 1, 2) = sign(ctm[2]);
      _os[0] = ctm[3];
      _os[1] = ctm[4];
      _os[2] = ctm[5];
   }

//--Store into compact storage

   void compact(int *const ctm) const
   {
      ctm[0] = _tm[0] + 2*_tm[3] + 3*_tm[6];
      ctm[1] = _tm[1] + 2*_tm[4] + 3*_tm[7];
      ctm[2] = _tm[2] + 2*_tm[5] + 3*_tm[8];
      ctm[3] = _os[0];
      ctm[4] = _os[1];
      ctm[5] = _os[2];
   }

//--Transform an index (matrix-vector multiplication)

   void transform(int *const index) const
   {
      const int i = index[0];
      const int j = index[1];
      const int k = index[2];
      index[0] = _tm[0]*i + _tm[1]*j + _tm[2]*k + _os[0];
      index[1] = _tm[3]*i + _tm[4]*j + _tm[5]*k + _os[1];
      index[2] = _tm[6]*i + _tm[7]*j + _tm[8]*k + _os[2];
   }

   void transform(int &i, int &j, int &k) const
   {
      const int ii = i;
      const int jj = j;
      const int kk = k;
      i = _tm[0]*ii + _tm[1]*jj + _tm[2]*kk + _os[0];
      j = _tm[3]*ii + _tm[4]*jj + _tm[5]*kk + _os[1];
      k = _tm[6]*ii + _tm[7]*jj + _tm[8]*kk + _os[2];
   }

//--Transpose

   void transpose()
   {
      const int tmp1 = _tm[1];
      _tm[1] = _tm[3];
      _tm[3] = tmp1;
      const int tmp2 = _tm[2];
      _tm[2] = _tm[6];
      _tm[6] = tmp2;
      const int tmp3 = _tm[5];
      _tm[5] = _tm[7];
      _tm[7] = tmp3;
   }

//--Reverse the transformation and offset to the other direction

   void reverse()
   {
      int ros[3];
      ros[0] = -_os[0];
      ros[1] = -_os[1];
      ros[2] = -_os[2];
      transpose();
      clear_offset();
      transform(ros);
      set_offset(ros);
   }

//--Clear the offset

   void clear_offset()
   {
      _os[0] = 0;
      _os[1] = 0;
      _os[2] = 0;
   }


/*==============================================================================
 * Private member functions
 *============================================================================*/

  private:

//--Constructors

   // From two transformation matrices.  This is used by the blkc library.
   CoordTransform(const CoordTransform &ct1, const CoordTransform &ct2)
   {
      matmul(ct1._tm, ct2._tm);
      clear_offset();
   }

//--Set the offset

   void set_offset(const int *const os)
   {
      _os[0] = os[0];
      _os[1] = os[1];
      _os[2] = os[2];
   }

//--Matrix multiplication

   void matmul(const int *const tm1, const int *const tm2)
   {
      _tm[0] = tm1[0]*tm2[0] + tm1[1]*tm2[3] + tm1[2]*tm2[6];
      _tm[1] = tm1[0]*tm2[1] + tm1[1]*tm2[4] + tm1[2]*tm2[7];
      _tm[2] = tm1[0]*tm2[2] + tm1[1]*tm2[5] + tm1[2]*tm2[8];
      _tm[3] = tm1[3]*tm2[0] + tm1[4]*tm2[3] + tm1[5]*tm2[6];
      _tm[4] = tm1[3]*tm2[1] + tm1[4]*tm2[4] + tm1[5]*tm2[7];
      _tm[5] = tm1[3]*tm2[2] + tm1[4]*tm2[5] + tm1[5]*tm2[8];
      _tm[6] = tm1[6]*tm2[0] + tm1[7]*tm2[3] + tm1[8]*tm2[6];
      _tm[7] = tm1[6]*tm2[1] + tm1[7]*tm2[4] + tm1[8]*tm2[7];
      _tm[8] = tm1[6]*tm2[2] + tm1[7]*tm2[5] + tm1[8]*tm2[8];
   }

//--Sign for transform

   int sign(const int x) const
   {
      return ( x < 0 ) ? -1 : 1;
   }


/*==============================================================================
 * Private operators
 *============================================================================*/

//--Index the matrix

   int &operator()(const int i, const int j) 
   {
      return _tm[3*i + j];
   }
   const int &operator()(const int i, const int j) const
   {
      return _tm[3*i + j];
   }

//--Assign constant to the matrix

   void operator=(const int val)
   {
      _tm[0] = val;
      _tm[1] = val;
      _tm[2] = val;
      _tm[3] = val;
      _tm[4] = val;
      _tm[5] = val;
      _tm[6] = val;
      _tm[7] = val;
      _tm[8] = val;
   }


/*==============================================================================
 * Friends
 *============================================================================*/

   friend class BlkC::BlockConnectivity;

};

#endif
