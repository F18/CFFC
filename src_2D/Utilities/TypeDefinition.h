/*!\file TypeDefinition.h
  \brief Header file defining different useful data types (e.g. template traits, runtime exception types etc.) */

#ifndef _TYPEDEFINITION_INCLUDED
#define _TYPEDEFINITION_INCLUDED

/* Include CFFC header files */
// None so far

/* Include required C++ libraries. */
#include <vector>
#include <stdexcept>

/* Using std namespace functions */
using std::vector;

/*!\var enum SpaceType {OneD=1, TwoD=2, ThreeD=3}
\brief Define the space dimensions.

One-dimension:    OneD <br>
Two-dimensions:   TwoD <br>
Three-dimensions: ThreeD 
*/
enum SpaceType {OneD=1, TwoD=2, ThreeD=3};

/* Type of polynomials */
enum PolynomOrder {CONSTANT = 0, LINEAR = 1, QUADRATIC = 2, CUBIC = 3, QUARTIC = 4, QUINTIC = 5 };

/* Index type --> vector for storing the indeces of the cells that are part of the stencil  */
typedef vector<int> IndexType;

// NumberOfVariables inside the solution class traits;
template<class T>
class SolutionParameters;

template<>
class SolutionParameters<double> {
 public:
  static const int NUM_OF_VARIABLES = 1; 
};



/***************** Type Definition *****************/
typedef double (* FunctionType1D) (const double);
typedef double (* FunctionType2D) (const double, const double);
typedef double (* FunctionType3D) (const double, const double, const double);
typedef double (* SuperFunctionType1D) (const double, const FunctionType1D, const FunctionType1D);

/***************** Exception Definitions *************/
/**
 * Exception to throw when the number of iteration
 * exceeds the maximum value
 */
struct maximum_exceeded : public std::logic_error
{
  maximum_exceeded() : std::logic_error("maximum number of iteration exceeded") {};
};

/**
 * Exception to throw when the required accuracy is
 * less than the machine accuracy
 */
struct machine_accuracy_exceeded : public std::logic_error
{
  machine_accuracy_exceeded(): std::logic_error("the values are beyond the accuracy of the machine") {};
};

/**
 * Exception to throw when the relative length of the interval for integration is
 * less than the machine accuracy
 */
struct too_short_interval : public std::logic_error
{
  too_short_interval(): std::logic_error("the dimension of the interval for integration is less than the machine accuracy") {};
};

/**
 * Exception to throw when a NULL pointer is not expected
 */

struct NullPointerException : public std::logic_error
{
 public:
   NullPointerException() : std::logic_error("Error: NULL pointer passed -> invalid argument"){ };
};


#endif
