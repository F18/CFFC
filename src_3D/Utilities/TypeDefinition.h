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
typedef enum {OneD=1, TwoD=2, ThreeD=3}  SpaceType;

/* Type of polynomials */
typedef enum {CONSTANT = 0, LINEAR = 1, QUADRATIC = 2, CUBIC = 3, QUARTIC = 4, QUINTIC = 5}  PolynomOrder;

/* Index type --> vector for storing the indexes of the cells that are part of the stencil  */
typedef vector<int> IndexType;

/* Short Index type --> vector for storing the indexes of type 'short int' */
typedef vector<short int> ShortIndexType;

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
 * The base for all CFFC exceptions.
 */
struct cffc_error : public std::exception
{
    cffc_error(const std::string& msg)
        : err_msg(msg)
    {
    }
    
    ~cffc_error() throw()
    {
    }
    
    const char* what() const throw()
    {
        return err_msg.c_str();
    }
    
private:

    std::string err_msg;
};

/**
 * Exception to be thrown when the number of iteration
 * exceeds the maximum allowed value
 */
struct MaximumIterationsExceeded : public cffc_error
{
  MaximumIterationsExceeded(const std::string& msg = "maximum number of iteration exceeded")
    : cffc_error(msg)
  {
  }

  ~MaximumIterationsExceeded() throw()
  {
  }
};

/**
 * Exception to be thrown when the required accuracy is
 * less than machine accuracy
 */
struct MachineAccuracyExceeded : public cffc_error
{
  MachineAccuracyExceeded(const std::string& msg = "machine accuracy exceeded")
    : cffc_error(msg)
  {
  }

  ~MachineAccuracyExceeded() throw()
  {
  }
};

/**
 * Exception to be thrown when the interval length 
 * (i.e. distance between 2 points) is less than 
 * machine accuracy
 */
struct TooShortInterval : public cffc_error
{
  TooShortInterval(const std::string& msg = "interval length shorter than machine accuracy")
    : cffc_error(msg)
  {
  }

  ~TooShortInterval() throw()
  {
  }
};

/**
 * Exception to be thrown when a not expected NULL is encountered
 */
struct ArgumentNullException : public cffc_error
{
  ArgumentNullException(const std::string& msg = "NULL pointer encountered => invalid argument")
    : cffc_error(msg)
  {
  }

  ~ArgumentNullException() throw()
  {
  }
};


#endif
