/*!\file HighOrder1D.h
  \brief Regression tests for functions prototyped in NumericalLibrary.h */

#ifndef _HIGHORDER_1D_INCLUDED
#define _HIGHORDER_1D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */

/* Include CFFC header files */
#include "TaylorDerivatives.h"
#include "../Math/Matrix.h"
#include "CENO_ExecutionMode.h"
#include "../Grid/Cell1D.h"


template<class SOLN_STATE>
class HighOrder1D{
  
public: 

  typedef SOLN_STATE Soln_State;

  /* Use TaylorDerivativesContainer_1D */
  typedef TaylorDerivativesContainer<OneD,Soln_State> DerivativesContainer;
  typedef TaylorDerivativesContainer<OneD,double> GeometricMoments;
  typedef typename DerivativesContainer::Derivative  Derivative;
  typedef Cell1D_Uniform GeometryType;
  
  // Constructors
  HighOrder1D(void);

  // Copy constructor
  HighOrder1D( const HighOrder1D & rhs);

  // Destructor
  ~HighOrder1D(void){ Geom = NULL;};

  // Assignment operator
  HighOrder1D<Soln_State> & operator=(const HighOrder1D<Soln_State> & rhs);

  /* Field access */
  const DerivativesContainer & CellDeriv(void) const {return TD;}
  DerivativesContainer & CellDeriv(void) {return TD;}
  const Soln_State & CellDeriv(int i) const {return TD(i);}
  Soln_State & CellDeriv(int i) {return TD(i);}
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

  const Soln_State & CellSmoothnessIndicator() const { return SI;}
  Soln_State & CellSmoothnessIndicator(){ return SI;}
  const double & CellSmoothnessIndicator(const int VarPosition)const{ return SI[VarPosition];}
  double & CellSmoothnessIndicator( const int VarPosition){ return SI[VarPosition];}

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

  // Cell geometry
  GeometryType* Geometry(void){return Geom;}
  const GeometryType* Geometry(void) const {return Geom;}
  GeometryType & CellGeometry(void){return *Geom;}
  const GeometryType & CellGeometry(void) const {return *Geom;}
  const double & CellCenter(void) const {return Geom->x;}
  const double & CellDelta (void) {return Geom->dx;}
  void SetGeometryPointer(GeometryType & Cell){ Geom = &Cell;}

  /* Operating functions */
  void SetRings(void);
  void ComputeGeometricCoefficients(void);
  double SolutionAtCoordinates(const double & X_Coord, const unsigned parameter){
    return TD.ComputeSolutionFor(X_Coord - CellCenter())[parameter];
  }
  Soln_State SolutionAtCoordinates(const double & X_Coord){
    return TD.ComputeSolutionFor(X_Coord - CellCenter());
  }

  /* Friend functions */
  friend bool operator== (const HighOrder1D<Soln_State> & left,
			  const HighOrder1D<Soln_State> & right);

  friend bool operator!= (const HighOrder1D<Soln_State> & left,
			  const HighOrder1D<Soln_State> & right);

  friend std::ostream & operator<< (std::ostream & os, const HighOrder1D<Soln_State> & Obj);
  friend std::istream & operator>> (std::istream & os, HighOrder1D<Soln_State> & Obj);


protected:
  
private:
  DerivativesContainer TD;  	// High-order TaylorDerivatives
  GeometricMoments GeomCoeff;   // The integrals of the geometric moments with respect to the centroid
  Soln_State SI;               // The values of the smoothness indicator calculated for each variable
  std::vector<short int> LimitedCell; // Monotonicity flag: Values --> OFF - high-order reconstruction
                                      //                                ON - limited linear reconstruction
  int rings;                    // Number of rings used to generate the reconstruction stencil
  int RecOrder;                 // Reconstruction Order

  /* Create storage for the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix CENO_LHS;
  ColumnVector CENO_Geometric_Weights;

  // Associate this reconstruction to a certain cell
  GeometryType* Geom;    // Pointer to cell geometry


  // Set Geometry_Ptr to point to the same geometry as the current Geom pointer
  void set_geometry_pointer(GeometryType* & Geometry_Ptr) const { Geometry_Ptr = Geom; }
};

// Constructors
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(void):TD(), GeomCoeff(), Geom(NULL),
					   CENO_LHS(), CENO_Geometric_Weights() {
  RecOrder = 0;
  rings = 0;
  LimitedCell.reserve(Soln_State::NumberOfVariables);
}

// Copy constructor 
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(const HighOrder1D<SOLN_STATE> & rhs): Geom(NULL)
{
  TD = rhs.CellDeriv();
  GeomCoeff = rhs.GeomCoeff;
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();
  RecOrder = rhs.CellRecOrder();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (rhs.GeomWeights().size() != 0)){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);
}

// Assignment operator
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE> & HighOrder1D<SOLN_STATE>::operator=(const HighOrder1D<SOLN_STATE> & rhs){

  // Handle self-assignment:
  if (this == & rhs) return *this;

  TD = rhs.CellDeriv();
  GeomCoeff = rhs.GeomCoeff;
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();
  RecOrder = rhs.CellRecOrder();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (rhs.GeomWeights().size() != 0)){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);

  return *this;
}



// Specializations
template<> inline
HighOrder1D<double>::HighOrder1D(void):TD(), GeomCoeff(), Geom(NULL),
				       CENO_LHS(), CENO_Geometric_Weights() {
  RecOrder = 0;
  rings = 0;
  LimitedCell.reserve(1);
}


#endif
