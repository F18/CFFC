/*!\file HighOrder1D.h
  \brief Regression tests for functions prototyped in NumericalLibrary.h */

#ifndef _HIGHORDER_1D_INCLUDED
#define _HIGHORDER_1D_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using std::ostream;
using std::istream;
using std::vector;

/* Include CFFC header files */
#include "TaylorDerivatives.h"
#include "../Math/Matrix.h"
#include "../Math/NumericalLibrary.h"
#include "CENO_ExecutionMode.h"
#include "../Grid/Cell1D.h"

/*********************************
 * Declare the HighOrder1D class *
 ********************************/
template <class SOLN_STATE> 
class HighOrder1D;

/************************************************
 *     Friend Functions : HighOrder1D           *
 ************************************************/
template<class SOLN_STATE>
bool operator== (const HighOrder1D<SOLN_STATE> & left,
		 const HighOrder1D<SOLN_STATE> & right);

template<class SOLN_STATE>
bool operator!= (const HighOrder1D<SOLN_STATE> & left,
		 const HighOrder1D<SOLN_STATE> & right);

template<class SOLN_STATE>
ostream & operator<< (ostream & os, const HighOrder1D<SOLN_STATE> & Obj);

template<class SOLN_STATE>
istream & operator>> (istream & os, HighOrder1D<SOLN_STATE> & Obj);


/*!
 * \class HighOrder1D
 *
 * \brief Template class for high-order variables in 1D.
 *********************************************************/
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
  HighOrder1D(int ReconstructionOrder, GeometryType & Cell);

  // Copy constructor
  HighOrder1D( const HighOrder1D & rhs);

  // Destructor
  ~HighOrder1D(void){ Geom = NULL;};

  // Assignment operator
  HighOrder1D<Soln_State> & operator=(const HighOrder1D<Soln_State> & rhs);

  /* Field access */
  const DerivativesContainer & CellDeriv(void) const {return TD;}
  DerivativesContainer & CellDeriv(void) {return TD;}
  const Soln_State & CellDeriv(const int & p1) const {return TD(p1);}
  Soln_State & CellDeriv(const int & p1) {return TD(p1);}
  const double & CellDeriv(const int & p1, const int & Variable) const { return TD(p1)[Variable];}
  double & CellDeriv(const int & p1, const int & Variable) { return TD(p1)[Variable];}
  Derivative & CellDeriv_InPosition(const int & position) {return TD(position,true,true,true);}
  const Derivative & CellDeriv_InPosition(const int & position) const {return TD(position,true,true,true);}

  const int NumberOfTaylorDerivatives() const {return TD.size();}

  const GeometricMoments & CellGeomCoeff() const { return GeomCoeff; }
  GeometricMoments & CellGeomCoeff() {return GeomCoeff;}
  const double & CellGeomCoeff(const int & p1) {return GeomCoeff(p1);}
  const double & CellGeomCoeff_InPosition(const int & position) {return GeomCoeff(position,true,true,true).D();}

  const int & CellRings() const {return rings;}
  const int & CellRecOrder() const {return TD.RecOrder();}

  /* Monotonicity variables --> high-order */
  const vector<short int> & CellInadequateFit() const { return LimitedCell;}
  vector<short int> & CellInadequateFit(){ return LimitedCell;}
  const short int & CellInadequateFit( const int VarPosition) const { return LimitedCell[VarPosition-1];}
  short int & CellInadequateFit( const int VarPosition){ return LimitedCell[VarPosition-1];}

  const Soln_State & CellSmoothnessIndicator() const { return SI;}
  Soln_State & CellSmoothnessIndicator(){ return SI;}
  const double & CellSmoothnessIndicator(const int & VarPosition)const{ return SI[VarPosition];}
  double & CellSmoothnessIndicator( const int & VarPosition){ return SI[VarPosition];}

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
  void InitializeMonotonicityFlag(void);

  double SolutionAtCoordinates(const double & X_Coord, const unsigned parameter){
    return TD.ComputeSolutionFor(X_Coord - CellCenter())[parameter];
  }
  Soln_State SolutionAtCoordinates(const double & X_Coord){
    return TD.ComputeSolutionFor(X_Coord - CellCenter());
  }
  /* Computation of right and left states (left and right are relative to the position of the centroid) */
  Soln_State left_state(void);
  Soln_State right_state(void);

  // Integrate over the domain of the geometry associated with this high-order solution
  template<typename FO, class ReturnType>
  ReturnType IntegrateOverTheCell(const FO FuncObj, const int & digits, ReturnType _dummy_param);

  /* Friend functions */
  friend bool operator== <Soln_State> (const HighOrder1D<Soln_State> & left,
				       const HighOrder1D<Soln_State> & right);

  friend bool operator!= <Soln_State> (const HighOrder1D<Soln_State> & left,
				       const HighOrder1D<Soln_State> & right);

  friend ostream & operator<< <Soln_State> (ostream & os, const HighOrder1D<Soln_State> & Obj);
  friend istream & operator>> <Soln_State> (istream & os, HighOrder1D<Soln_State> & Obj);


protected:
  
private:
  DerivativesContainer TD;  	// High-order TaylorDerivatives
  GeometricMoments GeomCoeff;   // The integrals of the geometric moments with respect to the centroid
  Soln_State SI;                // The values of the smoothness indicator calculated for each variable
  vector<short int> LimitedCell; // Monotonicity flag: Values --> OFF - high-order reconstruction
                                 //                               ON - limited linear reconstruction
  int rings;                    // Number of rings used to generate the reconstruction stencil

  /* Create storage for the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix CENO_LHS;
  ColumnVector CENO_Geometric_Weights;

  // Associate this reconstruction to a certain cell
  GeometryType* Geom;    // Pointer to cell geometry


  // Set Geometry_Ptr to point to the same geometry as the current Geom pointer
  void set_geometry_pointer(GeometryType* & Geometry_Ptr) const { Geometry_Ptr = Geom; }
};

/****************************************************
 * Implement the HighOrder1D class member functions *
 ***************************************************/

//! Default Constructor
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(void):TD(0), GeomCoeff(0), SI(0),
					   CENO_LHS(), CENO_Geometric_Weights(), Geom(NULL) {
  // set rings
  SetRings();

  // initialize monotonicity flag
  InitializeMonotonicityFlag();

  // set the value of the geometric moment
  GeomCoeff(0) = 1.0;
}

//! Main Constructor
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(int ReconstructionOrder, GeometryType & Cell):
  TD(ReconstructionOrder), GeomCoeff(ReconstructionOrder),
  SI(), CENO_LHS(), CENO_Geometric_Weights(){

  // set geometry pointer
  SetGeometryPointer(Cell);

  // set rings
  SetRings();
  
  // initialize monotonicity flag
  InitializeMonotonicityFlag();

  // compute geometric moments for the provided geometry
  ComputeGeometricCoefficients();
}

//! Copy constructor 
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(const HighOrder1D<SOLN_STATE> & rhs): Geom(NULL)
{
  TD = rhs.CellDeriv();
  GeomCoeff = rhs.GeomCoeff;
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (rhs.GeomWeights().size() != 0)){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);
}

//! Assignment operator
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE> & HighOrder1D<SOLN_STATE>::operator=(const HighOrder1D<SOLN_STATE> & rhs){

  // Handle self-assignment:
  if (this == & rhs) return *this;

  TD = rhs.CellDeriv();
  GeomCoeff = rhs.GeomCoeff;
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (!rhs.GeomWeights().null()) ){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);

  return *this;
}

//! SetRings()
/*! Set the number of rings around the current cell
 * which will be used to form the supporting stencil
 * of the reconstruction. 
 * In 1D, each ring adds two neighbours to the stencil.
 * TD.size() --> gives the number of TaylorDerivatives,
 * which is a function of the order of reconstruction
 */
template<class SOLN_STATE> inline
void HighOrder1D<SOLN_STATE>::SetRings(void){
 
  switch(TD.size()){
  case 1: rings = 0; break;	// piecewise constant
  case 2: rings = 1; break;	// piecewise linear
  case 3: rings = 2; break;	// piecewise quadratic
  case 4: rings = 2; break;	// piecewise cubic
  case 5: rings = 3; break;	// piecewise quartic
  case 6: rings = 3; break;	// piecewise quintic
  case 7: rings = 4; break;	// piecewise 6th-order
  case 8: rings = 4; break;	// piecewise 7th-order
  case 9: rings = 5; break;	// piecewise 8th-order
  default: rings = -1;
  }
}

//! ComputeGeometricCoefficients()
/*! 
 * 
 * Computes the geometric coefficients by computing integral      
 * over the volume of the cell of the polynomial function having  
 * the form                                                       
 *       ((x-xCC)^n) and divide by length dx                       
 */
template<class SOLN_STATE>
void HighOrder1D<SOLN_STATE>::ComputeGeometricCoefficients(void){

  /* Create the polynomial function for the cell */
  GeneralizedPolynomialFunctionOfOneVariable Polynom(0,CellCenter());
  GeomCoeff(0) = 1.0;
  for (int p1=1; p1<=GeomCoeff.RecOrder(); ++p1){
    Polynom.ChangePowerTo(p1);
    GeomCoeff(p1) = IntegrateOverTheCell(Polynom,14,GeomCoeff(p1))/CellDelta();
  }
}

//! Initialize monotonicity flag
template<class SOLN_STATE> inline
void HighOrder1D<SOLN_STATE>::InitializeMonotonicityFlag(void){
  LimitedCell.reserve(Soln_State::NumberOfVariables);
  for (int i = 0; i <= Soln_State::NumberOfVariables - 1; ++i){
    LimitedCell[i] = OFF; /* initialize the flags to OFF (smooth solution)*/
  }
}

//! Integrate over the cell
template<class SOLN_STATE> template<typename FO, class ReturnType> inline
ReturnType HighOrder1D<SOLN_STATE>::IntegrateOverTheCell(const FO FuncObj,
							 const int & digits,
							 ReturnType _dummy_param){
  return AdaptiveGaussianQuadrature(FuncObj,
				    CellCenter() - 0.5* CellDelta(),
				    CellCenter() + 0.5* CellDelta(),
				    _dummy_param,digits);
}


//! Compute the solution at the left cell interface
template<class SOLN_STATE> inline
SOLN_STATE HighOrder1D<SOLN_STATE>::left_state(void){
  return SolutionAtCoordinates(CellCenter() - 0.5* CellDelta());
}

//! Compute the solution at the right cell interface
template<class SOLN_STATE> inline
SOLN_STATE HighOrder1D<SOLN_STATE>::right_state(void){
  return SolutionAtCoordinates(CellCenter() + 0.5* CellDelta());
}

// Friend functions
//! operator== 
template<class SOLN_STATE> inline
bool operator== (const HighOrder1D<SOLN_STATE> & left, const HighOrder1D<SOLN_STATE> & right){
   
  if ( CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (!(left.GeomWeights().null() && right.GeomWeights().null())) ){
    return ( (left.CellDeriv() == right.CellDeriv()) && 
	     (left.CellGeomCoeff() == right.CellGeomCoeff()) &&
	     (left.CellSmoothnessIndicator() == right.CellSmoothnessIndicator()) &&
	     (left.CellInadequateFit() == right.CellInadequateFit()) &&
	     (left.CellRings() == right.CellRings()) && 
	     (left.CellGeometry() == right.CellGeometry()) && 
	     (left.LHS() == right.LHS()) &&
	     (left.GeomWeights() == right.GeomWeights()) 
	     );
  } else {
    // There is no memory allocated for pseudo-inverse data in both objects
    return ( (left.CellDeriv() == right.CellDeriv()) && 
	     (left.CellGeomCoeff() == right.CellGeomCoeff()) &&
	     (left.CellSmoothnessIndicator() == right.CellSmoothnessIndicator()) &&
	     (left.CellInadequateFit() == right.CellInadequateFit()) &&
	     (left.CellRings() == right.CellRings()) && 
	     (left.CellGeometry() == right.CellGeometry()) );
  }
}

//! operator!= 
template<class SOLN_STATE> inline
bool operator!= (const HighOrder1D<SOLN_STATE> & left, const HighOrder1D<SOLN_STATE> & right){
  return !(left == right);
}

//! operator<<
template<class SOLN_STATE> inline
ostream & operator<< (ostream & os, const HighOrder1D<SOLN_STATE> & Obj){

  os.setf(ios::skipws,ios::scientific);
  os << Obj.CellDeriv() 
     << Obj.CellGeomCoeff();
  os.width(4);
  os << Obj.CellSmoothnessIndicator() << endl;
  os.width(4);
  os << Obj.CellRings();
  os.unsetf(ios::skipws);
  os.unsetf(ios::scientific);
  return os;
}

//! operator>>
template<class SOLN_STATE> inline
istream & operator>> (istream & os, HighOrder1D<SOLN_STATE> & Obj){

  os.setf(ios::skipws);
  os >> Obj.CellDeriv()
     >> Obj.CellGeomCoeff()
     >> Obj.CellSmoothnessIndicator();

  // Set rings
  Obj.SetRings();
  os.unsetf(ios::skipws);
  return os;
}


// Specializations
template<> inline
void HighOrder1D<double>::InitializeMonotonicityFlag(void){
  LimitedCell.reserve(1);
  LimitedCell[0] = OFF;
}

template<> inline
double & HighOrder1D<double>::CellDeriv(const int & p1, const int & Variable) {
 return TD(p1);
}

template<> inline
const double & HighOrder1D<double>::CellDeriv(const int & p1, const int & Variable) const {
 return TD(p1);
}

template<> inline
const double & HighOrder1D<double>::CellSmoothnessIndicator(const int & VarPosition) const {
  return SI;
}

template<> inline
double & HighOrder1D<double>::CellSmoothnessIndicator(const int & VarPosition) {
  return SI;
}

template <> inline
double HighOrder1D<double>::SolutionAtCoordinates(const double & X_Coord, const unsigned parameter){
  return TD.ComputeSolutionFor(X_Coord - CellCenter());
}


#endif
