/*!\file RandomGenerator.h
  \brief Header file defining useful numerical methods/objects related to random number generators. */

#ifndef _RANDOM_NUMBER_GENERATOR_INCLUDED
#define _RANDOM_NUMBER_GENERATOR_INCLUDED

/* Include required C++ libraries. */
#include <cmath>

/* Using std namespace functions */
// None

/* Include CFFC header files */
#include "../Utilities/Utilities.h"

/*!
 * \class RandomGen
 * 
 * Implementation of the highest quality recommended 
 * generator in Numerical Recipes Third Ed. (pp. 342)
 * Once a generator is initialized it will maintain the
 * internal state between calls.
 * If you need a random integer between 1 and n (inclusive),
 * 
 */
class RandomGen {
  
public:

  typedef unsigned long long int Ullong;

  //! @brief Constructor
  RandomGen( Ullong seed);

  //! Return 64-bit random integer
   Ullong int64(void);

  //! Return random double-precision floating value in the range of 0. to 1.
  double doub(void) { return 5.42101086242752217E-20 * int64(); }

  //! Return 32-bit random integer
  unsigned int int32(void) { return (unsigned int)int64(); }

private:
   Ullong u, v, w;
  
  // Default private constructor
  RandomGen(void);
};

/*!
 * Constructor with seed.
 *
 * \note Call it with any integer seed exept 4101842887655102017LL!
 */
inline RandomGen::RandomGen(Ullong seed):
  v(4101842887655102017LL), w(1) 
{
  u = seed ^ v;  int64();
  v = u;  int64();
  w = v;  int64();
}


/*!
 * \class Expondev
 * 
 * Generator for exponential deviates (Numerical Recipes Third Ed. (pp. 362) )
 * 
 */
class Expondev : RandomGen {

public:

  Expondev(Ullong seed) : RandomGen(seed), beta(1.0) {};
  Expondev(double bbeta, Ullong seed) : RandomGen(seed), beta(bbeta) {};
  
  double dev(void) {
    double u;
    do u = doub(); while (u == 0.);
    return -log(u)/beta;
  }

private:
  double beta;
};

#endif
