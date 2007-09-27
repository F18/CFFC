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

  // Constructors
  HighOrder1DState(void){}

  // Copy constructor



  // Destructor
  ~HighOrder1DState(void){}

  /* Field access */
  const DerivativesContainer & CellDeriv(void) const {return TD;}
  DerivativesContainer & CellDeriv(void) {return TD;}
  const Soln_pState & CellDeriv(int i) const {return TD(i);}
  Soln_pState & CellDeriv(int i) {return TD(i);}
  const double & CellDeriv(const int i, const int Variable) const { return TD(i)[Variable];}
  double & CellDeriv(const int i, const int Variable) { return TD(i)[Variable];}
  Derivative & CellDeriv(const int position, const bool, const bool) {return TD(position,true,true,true);}
  const Derivative & CellDeriv(const int position, const bool, const bool) const {return TD(position,true,true,true);}

  const int NumberOfTaylorDerivatives() const {return TD.size();}

  const GeometricMoments & CellGeomCoeff() const { return GeomCoeff; }
  GeometricMoments & CellGeomCoeff() {return GeomCoeff;}
  const double & CellGeomCoeff(const int & p1) {return GeomCoeff(p1);}
  const double & CellGeomCoeff(const int position, const bool) {return GeomCoeff(position,true,true,true).D();}

  const int & CellRings() const {return rings;}
  const int & CellRecOrder() const {return RecOrder;}

  /* Monotonicity variables --> high-order */
  const std::vector<short int> & CellInadequateFit() const { return LimitedCell;}
  std::vector<short int> & CellInadequateFit(){ return LimitedCell;}
  const short int & CellInadequateFit( const int VarPosition) const { return LimitedCell[VarPosition-1];}
  short int & CellInadequateFit( const int VarPosition){ return LimitedCell[VarPosition-1];}

  const Soln_pState & CellSmoothnessIndicator() const { return SI;}
  Soln_pState & CellSmoothnessIndicator(){ return SI;}
  const double & CellSmoothnessIndicator(const int VarPosition)const{ return SI[VarPosition];}
  double & CellSmoothnessIndicator( const int VarPosition){ return SI[VarPosition];}

#ifdef _CENO_SPEED_EFFICIENT
  /* Access the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix & LHS(void) {return CENO_LHS;}
  const DenseMatrix & LHS(void) const {return CENO_LHS;}
  double & LHS(const int & IndexI, const int & IndexJ) {return CENO_LHS(IndexI,IndexJ);}
  const double & LHS(const int & IndexI, const int & IndexJ) const {return CENO_LHS(IndexI,IndexJ);}
  /* Access the geometric weights assigned to the cells in the stencil */
  ColumnVector & GeomWeights(void){return CENO_Geometric_Weights;}
  const ColumnVector & GeomWeights(void) const {return CENO_Geometric_Weights;}
  double & GeomWeights(const int & Index){return CENO_Geometric_Weights(Index);}
  const double & GeomWeights(const int & Index) const {return CENO_Geometric_Weights(Index);}
#endif //_CENO_SPEED_EFFICIENT

  /* Operating functions */
  void SetRings(void);
  void ComputeGeometricCoefficients(void);
  double SolutionAtCoordinates(const double & X_Coord, const double & CellCenter, const unsigned parameter){
    return TD.ComputeSolutionFor(X_Coord - CellCenter)[parameter];
  }

protected:
  
  
private:
  DerivativesContainer TD;  	// High-order TaylorDerivatives
  GeometricMoments GeomCoeff;   // The integrals of the geometric moments with respect to the centroid
  Soln_pState SI;               // The values of the smoothness indicator calculated for each variable
  std::vector<short int> LimitedCell; // Monotonicity flag: Values --> OFF - high-order reconstruction
                                      //                                ON - limited linear reconstruction
  int rings;                    // Number of rings used to generate the reconstruction stencil
  int RecOrder;                 // Reconstruction Order

#ifdef _CENO_SPEED_EFFICIENT
  /* Create storage for the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix CENO_LHS;
  ColumnVector CENO_Geometric_Weights;
#endif //_CENO_SPEED_EFFICIENT

};
