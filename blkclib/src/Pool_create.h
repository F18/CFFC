#ifndef _POOL_CREATE_INCLUDED
#define _POOL_CREATE_INCLUDED

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


/*******************************************************************************
 *
 * A quick pool for only creating objects.  Deletion is not supported.
 *
 ******************************************************************************/

namespace Pool
{

//--Class block holds a chunk of memory for type T

template <typename T>
class Block
{
  public:
   T *array;
   Block *prev;

   Block(const size_t size, Block *const _prev)
      :
      prev(_prev)
   {
      array = new T[size];
   }

   ~Block()
   {
      delete[] array;
   }

  private:
   Block(const Block&);
   Block &operator=(const Block&);
};

//--Class Create has interface for getting objects T

template <typename T>
class Create
{
  public:
   Create(const size_t _size)
      :
      size(_size),
      index(_size),
      block(0)
   { }

   ~Create() {
      while ( block ) {
         Block<T> *const prev = block->prev;
         delete block;
         block = prev;
      }
   }

   T *const get()
   {
      if ( index == size ) get_new_block();
      return &(block->array[index++]);
   }

  private:
   const size_t size;
   size_t index;
   Block<T> *block;

   void get_new_block() {
      Block<T> *const prev = block;
      block = new Block<T>(size, prev);
      index = 0;
   }

   Create(const Create&);
   Create &operator=(const Create&);
};

}  // namespace Pool

#endif
