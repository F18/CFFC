/*!\file HighOrder1DState.h
  \brief Regression tests for functions prototyped in NumericalLibrary.h */

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TaylorDerivatives.h"
#include "../Math/Matrix.h"


template<class SOLN_pSTATE, class SOLN_cSTATE>
class HighOrder1DState{
  
public: 

  typedef SOLN_pSTATE Soln_pState;
  typedef SOLN_cSTATE Soln_cState;

  /* Use TaylorDerivativesContainer_1D */
  typedef TaylorDerivativesContainer<OneD,Soln_pState> DerivativesContainer;
  typedef TaylorDerivativesContainer<OneD,double> GeometricMoments;
  typedef typename DerivativesContainer::Derivative  Derivative;

protected:
  
  
private:
  DerivativesContainer TD;  	// High-order TaylorDerivatives
  GeometricMoments GeomCoeff;   // The integrals of the geometric moments with respect to the centroid
  Soln_pState CharactVar;       // The characteristic variables
  Soln_pState SI;               // The values of the smoothness indicator calculated for each variable
  std::vector<short int> LimitedCell; // Monotonicity flag: Values --> OFF - high-order reconstruction
                                      //                                ON - limited linear reconstruction
#ifdef _CENO_SPEED_EFFICIENT
  /* Create storage for the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix CENO_LHS;
  ColumnVector CENO_Geometric_Weights;
#endif //_CENO_SPEED_EFFICIENT

};
