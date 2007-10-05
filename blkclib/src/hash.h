#ifndef _HASH_INCLUDED
#define _HASH_INCLUDED

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

//--FNV hashing parameters

#if defined(HAVE_64BIT_SIZE_T)
#define FNV_PRIME        1099511628211UL
#define FNV_OFFSET_BASIS 14695981039346656037UL
#else
#define FNV_PRIME        16777619UL
#define FNV_OFFSET_BASIS 2166136261UL
#endif

//--Hash FNV1a implemented via for loop.  "key" has size "len" bytes.

inline size_t hash_FNV1a(const void *const key, const int len)
{
  const unsigned char *p = static_cast<const unsigned char*>(key);
  size_t hash = FNV_OFFSET_BASIS;
  for( int n = len; n--; ) hash = (hash^static_cast<size_t>(*p++))*FNV_PRIME;
  return hash;
}

//--Hash FNV1a implemented via template-metaprogramming loop.  This should be
//--used if the length N is known at compile time.  "key" has size "N" bytes.
//--Use the entry point HashFNV1a<N>::eval(key).

template <int N> struct Hash1FNV1a {
  static size_t eval(size_t hash, const unsigned char *p)
  {
    return Hash1FNV1a<N-1>::eval((hash^static_cast<size_t>(*p))*FNV_PRIME, p+1);
  }
};

template <> struct Hash1FNV1a<1> {
  static size_t eval(size_t hash, const unsigned char *p)
  {
    return (hash^static_cast<size_t>(*p))*FNV_PRIME;
  }
};

// Entry point
template <int N> struct HashFNV1a {
  static size_t eval(const void *const key) 
  {
    size_t hash = FNV_OFFSET_BASIS;
    return Hash1FNV1a<N>::eval(hash, static_cast<const unsigned char*>(key));
  }
};

#endif
