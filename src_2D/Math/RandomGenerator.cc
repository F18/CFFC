#include "RandomGenerator.h"

/*!
 * Return 64-bit random integer.
 */
RandomGen::Ullong RandomGen::int64(void){
  
  u = ( u * 2862933555777941757LL + 
	7046029225438635308LL );
  v ^= v >> 17; 
  v ^= v << 31;
  v ^= v >> 8;
  w = 4294957665U * (w & 0xffffffff) + (w >> 32);

  Ullong x = u ^ (u << 21);
  x ^= x >> 35; 
  x ^= x << 4;

  return (x + v) ^ w;
}
