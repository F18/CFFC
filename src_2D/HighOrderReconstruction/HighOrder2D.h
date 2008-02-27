/*!\file HighOrder2D.h
  \brief Header file implementing the templated HighOrder2D class.*/

#ifndef _HIGHORDER_2D_INCLUDED
#define _HIGHORDER_2D_INCLUDED

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
#include "../Grid/HO_Cell2D.h"
//#include "ReconstructionSolvers2D.h"

/*********************************
 * Declare the HighOrder2D class *
 ********************************/
template <class SOLN_STATE> 
class HighOrder2D;

/************************************************
 *     Friend Functions : HighOrder2D           *
 ************************************************/
template<class SOLN_STATE>
bool operator== (const HighOrder2D<SOLN_STATE> & left,
		 const HighOrder2D<SOLN_STATE> & right);

template<class SOLN_STATE>
bool operator!= (const HighOrder2D<SOLN_STATE> & left,
		 const HighOrder2D<SOLN_STATE> & right);

template<class SOLN_STATE>
ostream & operator<< (ostream & os, const HighOrder2D<SOLN_STATE> & Obj);

template<class SOLN_STATE>
istream & operator>> (istream & os, HighOrder2D<SOLN_STATE> & Obj);


/*!
 * \class HighOrder2D
 *
 * \brief Template class for high-order variables in 2D.
 * \nosubgrouping
 *********************************************************/
template<class SOLN_STATE>
class HighOrder2D{
  
public: 

  //! @name Defined public types:
  //@{
  typedef SOLN_STATE Soln_State;
  typedef HighOrder2D<Soln_State> ClassType;
  typedef Grid2D_Quad_Block_HO GeometryType;
  typedef TaylorDerivativesContainer<TwoD,Soln_State> DerivativesContainer;
  typedef typename DerivativesContainer::Derivative  Derivative;
  typedef typename GeometryType::GeometricMoments GeometricMoments;
  typedef typename GeometryType::GeomMoment       GeomMoment;
  //! Type for monotonicity information of each reconstructed variable.
  typedef std::vector<short int> FlagType;
  //! Type for smoothness indicator associated with each reconstructed variable.
  typedef std::vector<double> SmoothnessIndicatorType;

  typedef double (ClassType::*MemberFunction_OneArgument_Type)(const double &);
  typedef double (ClassType::*MemberFunction_TwoArguments_Type)(const double &, const double &);
  typedef double (ClassType::*MemberFunction_TwoArguments_OneParameter_Type)(const double &, const unsigned );
  //@}
  
  //! @name Constructors:
  //@{
  HighOrder2D(void);		//!< Simple constructor 
  HighOrder2D(int ReconstructionOrder, GeometryType & Block); //!< Advanced constructor
  HighOrder2D( const HighOrder2D & rhs); //!< Copy constructor
  //@}
  

  //! @name Destructors:
  //@{
  ~HighOrder2D(void){ Geom = NULL;}
  void deallocate(void);
  //@} 

  HighOrder2D<Soln_State> & operator=(const HighOrder2D<Soln_State> & rhs); //!< Assignment operator

  //! @name Taylor derivatives:
  //@{
  //! Get the pointer to the double array of Taylor derivatives.
  DerivativesContainer ** const TaylorDeriv(void) const {return TD;}
  //! Get the pointer to the double array of Taylor derivatives.
  DerivativesContainer ** TaylorDeriv(void) {return TD;}
  //! Get the container of Taylor derivatives for cell (ii,jj) of the 2D block.
  const DerivativesContainer & CellTaylorDeriv(const int &ii, const int &jj) const {return TD[ii][jj];}
  //! Get the container of Taylor derivatives for cell (ii,jj) of the 2D block.
  DerivativesContainer & CellTaylorDeriv(const int &ii, const int &jj) {return TD[ii][jj];}
  //! Get the solution state of the Taylor derivative of cell (ii,jj) for the (p1,p2) powers.
  const Soln_State & CellTaylorDerivState(const int &ii, const int &jj,
					  const int & p1, const int & p2) const {return TD[ii][jj](p1,p2);}
  //! Get the solution state of the Taylor derivative of cell (ii,jj) for the (p1,p2) powers.
  Soln_State & CellTaylorDerivState(const int &ii, const int &jj,
				    const int & p1, const int & p2) {return TD[ii][jj](p1,p2);}
  //! Get the value of Taylor derivative of cell (ii,jj) for the (p1,p2) powers and the specified 'Variable'.
  const double & CellTaylorDerivValue(const int & ii, const int & jj, 
				      const int & p1, const int & p2, const int & Variable) const {
    return TD[ii][jj](p1,p2)[Variable];
  }
  //! Get the value of Taylor derivative of cell (ii,jj) for the (p1,p2) powers and the specified Variable.
  double & CellTaylorDerivValue(const int & ii, const int & jj,
				const int & p1, const int & p2, const int & Variable) { return TD[ii][jj](p1,p2)[Variable];}
  //! Get the Taylor derivative of cell (ii,jj) which is stored in the 'position' place (i.e. powers and value).
  Derivative & CellTaylorDeriv(const int & ii, const int & jj,
			       const int & position) {return TD[ii][jj](position);}
  //! Get the Taylor derivative of cell (ii,jj) which is stored in the 'position' place (i.e. powers and value).
  const Derivative & CellTaylorDeriv(const int & ii, const int & jj,
				     const int & position) const {return TD[ii][jj](position);}
  //! Get the number of Taylor derivatives in each cell container
  const int NumberOfTaylorDerivatives(void) const {return TD[0][0].size();}
  //@} 

  //! @name Geometric moments
  //@{
  //! Get the geometric coefficients of cell (ii,jj).
  const GeometricMoments & CellGeomCoeff(const int & ii, const int & jj) const { return Geom->CellGeomCoeff(ii,jj); }
  //! Get the geometric coefficients of cell (ii,jj) which is stored in the 'position' place.
  const GeomMoment & CellGeomCoeff(const int & ii, const int & jj,
				   const int & position) {return GeomCoeff(ii,jj,position);}
  //! Get the value of cell (ii,jj) geometric coefficient with x-power 'p1' and y-power 'p2'.
  const double & CellGeomCoeffValue(const int & ii, const int & jj,
				    const int & p1, const int & p2) {return Geom->GeomCoeffValue(ii,jj,p1,p2);}
  //! Get the cell (ii,jj) geometric coefficient which is store in the 'position' place (i.e. powers and value)
  const double & CellGeomCoeffValue(const int & ii, const int & jj,
				    const int & position) {return Geom->GeomCoeffValue(ii,jj,position);}
  //@}

  //! @name Reconstruction order and number of rings
  //@{
  //! Get the number of rings of neighbour cells around the reconstructed cell.
  const int & Rings(void) const {return rings;}
  //! Get the order of reconstruction set for this object.
  const int RecOrder(void) const {return TD[0][0].RecOrder();}
  //@}

  //! @name Monotonicity info for high-order
  //@{
  //! Get the pointer to the double array of monotonicity flags.
  FlagType ** const InadequateFit(void) const { return LimitedCell;}
  //! Get the monotonicity flags for reconstruction of cell (ii,jj)
  const FlagType & CellInadequateFit(const int & ii, const int & jj) const { return LimitedCell[ii][jj];}
  //! Get the monotonicity flags for reconstruction of cell (ii,jj)
  FlagType & CellInadequateFit(const int & ii, const int & jj){ return LimitedCell[ii][jj];}
  //! Get the monotonicity flag for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'. 
  const short int & CellInadequateFitValue(const int & ii, const int & jj,
					   const int VarPosition) const { return LimitedCell[ii][jj][VarPosition-1];}
  //! Get the monotonicity flag for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  short int & CellInadequateFitValue(const int & ii, const int & jj,
				     const int VarPosition){ return LimitedCell[ii][jj][VarPosition-1];}
  //@}

  //! @name Smoothness indicator
  //@{
  //! Get the pointer to the double array of reconstruction smoothness indicators.
  SmoothnessIndicatorType ** const SmoothnessIndicator(void) const { return SI;}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj)
  const SmoothnessIndicatorType & CellSmoothnessIndicator(const int & ii, const int & jj) const { return SI[ii][jj];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) 
  SmoothnessIndicatorType & CellSmoothnessIndicator(const int & ii, const int & jj){ return SI[ii][jj];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  const double & CellSmoothnessIndicatorValue(const int & ii, const int & jj,
					      const int & VarPosition) const { return SI[ii][jj][VarPosition-1];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  double & CellSmoothnessIndicatorValue(const int & ii, const int & jj,
					const int & VarPosition){ return SI[ii][jj]a[VarPosition-1];}
  //@}

  //! @name Pseudo-inverse of the LHS term in the CENO reconstruction
  //@{
  //! Get the pointer to the double array of reconstruction pseudo-inverse matrices.
  DenseMatrix ** const LHS_Inv(void) {return CENO_LHS;}
  //! Get the pseudo-inverse matrix for the reconstruction of cell (ii,jj)
  DenseMatrix & Cell_LHS_Inv(const int & ii, const int & jj) {return CENO_LHS[ii][jj];}
  //! Get the pseudo-inverse matrix for the reconstruction of cell (ii,jj)
  const DenseMatrix & Cell_LHS_Inv(const int & ii, const int & jj) const {return CENO_LHS[ii][jj];}
  //! Get the entry (IndexI,IndexJ) in the pseudo-inverse matrix for the reconstruction of cell (ii,jj)
  double & Cell_LHS_Inv_Value(const int & ii, const int & jj,
			      const int & IndexI, const int & IndexJ) {return CENO_LHS[ii][jj](IndexI,IndexJ);}
  //! Get the entry (IndexI,IndexJ) in the pseudo-inverse matrix for the reconstruction of cell (ii,jj)
  const double & Cell_LHS_Inv_Value(const int & ii, const int & jj,
				    const int & IndexI, const int & IndexJ) const {return CENO_LHS(IndexI,IndexJ);}
  //! Return true if the pseudo-inverse has been already computed, otherwise false.
  bool IsPseudoInversePreComputed(void) const { return (PseudoInverseFlag==ON) ? true:false; }
  //! Require update of the pseudo-inverse
  void MustUpdatePseudoInverse(void) {PseudoInverseFlag = OFF;}
  //@}

  //! @name Geometric weights assigned to the cells that are part of the stencil
  //@{
  //! Get the pointer to the double array of geometric weights used for k-exact reconstruction.
  ColumnVector ** const GeomWeights(void){return CENO_Geometric_Weights;}

  ColumnVector & GeomWeights(void){return CENO_Geometric_Weights;}
  const ColumnVector & GeomWeights(void) const {return CENO_Geometric_Weights;}
  double & GeomWeights(const int & Index){return CENO_Geometric_Weights(Index);}
  const double & GeomWeights(const int & Index) const {return CENO_Geometric_Weights(Index);}
  //@}

  //! @name Cell geometry
  //@{
  GeometryType* Geometry(void){return Geom;}
  const GeometryType* Geometry(void) const {return Geom;}
  GeometryType & CellGeometry(void){return *Geom;}
  const GeometryType & CellGeometry(void) const {return *Geom;}
  const Vector2D & CellCenter(void) const {return Geom->Xc;}
  const double & XCellCenter(void) const {return Geom->Xc.x;}
  const double & YCellCenter(void) const {return Geom->Xc.y;}
  void SetGeometryPointer(GeometryType & Block){ Geom = &Block;}
  void AssociateGeometry(GeometryType & Block);
  //@}

  //! @name Initialize container functions.
  //@{
  static int Nghost(int ReconstructionOrder); //!< return the required number of ghost cells
  void SetRings(void);
  int StencilSize(void) const;
  void ComputeGeometricCoefficients(void);
  void InitializeMonotonicityFlag(void);
  void InitializeSmoothnessIndicator(void);

  void InitializeVariable(int ReconstructionOrder, int ReconstructionMethod);
  template<typename InputParametersType>
  void InitializeVariable(const InputParametersType & IP);
  //@}

  //! @name Evaluate polynomial interpolant functions.
  //@{
  double SolutionAtCoordinates(const double & X_Coord, const double & Y_Coord, const unsigned & parameter){
    return TD.ComputeSolutionFor(X_Coord - XCellCenter(), Y_Coord - YCellCenter())[parameter];
  }
  Soln_State SolutionAtCoordinates(const double & X_Coord, const double & Y_Coord){
    return TD.ComputeSolutionFor(X_Coord - XCellCenter(), Y_Coord - YCellCenter());
  }
  //@}

  /*! @brief Integrate over the domain of the geometry associated with this high-order solution  */
  template<typename FO, class ReturnType>
  ReturnType IntegrateOverTheCell(const FO FuncObj, const int & digits, ReturnType _dummy_param);
  
  //! @name Reconstructions:
  //@{
  /*! @brief Compute the unlimited high-order solution reconstruction.  */
  template<class Soln_Block_Type>
  void ComputeUnlimitedSolutionReconstruction(Soln_Block_Type *SolnBlk,
					      const int &iCell,
					      const int &jCell,
					      const int ReconstructionMethod,
					      HighOrder2D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
					      &Soln_Block_Type::CellHighOrder);

  /*! @brief Compute the pseudo-inverse corresponding to the unlimited high-order solution reconstruction.  */
  template<class Soln_Block_Type>
  void ComputeReconstructionPseudoInverse(Soln_Block_Type *SolnBlk,
					  const int &iCell,
					  const int &jCell,
					  HighOrder2D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
					  &Soln_Block_Type::CellHighOrder);

  /*! @brief Compute the second (low-order) reconstruction in the CENO algorithm.  */
  template<class Soln_Block_Type>
  void ComputeLowOrderReconstruction(Soln_Block_Type *SolnBlk,
				     const int &iCell,
				     const int &jCell,
				     const int &Limiter);

  //@}

  //! @name CENO Analysis:
  //@{
  template<class Soln_Block_Type>
  void ComputeSmoothnessIndicator(Soln_Block_Type *SolnBlk,
				  const int &iCell,
				  const int &jCell,
				  HighOrder2D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
				  &Soln_Block_Type::CellHighOrder);
  //@}


  //! @name Error Evaluation:
  //@{
  /*! @brief Compute L1 norm of the solution error */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL1(Function_Object_Type FuncObj, const unsigned &parameter);

  double ComputeSolutionErrorL1(HighOrder2D<Soln_State> & Obj, const unsigned &parameter);

  /*! @brief Compute the L2 norm of the solution error */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL2(Function_Object_Type FuncObj, const unsigned &parameter);

  double ComputeSolutionErrorL2(const HighOrder2D<Soln_State> & Obj);
  //@}

  /* Friend functions */
  friend bool operator== <Soln_State> (const HighOrder2D<Soln_State> & left,
				       const HighOrder2D<Soln_State> & right);

  friend bool operator!= <Soln_State> (const HighOrder2D<Soln_State> & left,
				       const HighOrder2D<Soln_State> & right);

  friend ostream & operator<< <Soln_State> (ostream & os, const HighOrder2D<Soln_State> & Obj);
  friend istream & operator>> <Soln_State> (istream & os, HighOrder2D<Soln_State> & Obj);


protected:
  
private:
  DerivativesContainer **TD;  	   //!< High-order TaylorDerivatives
  vector<double> **SI;             //!< The values of the smoothness indicator calculated for each reconstructed variable
  vector<short int> **LimitedCell; //!< Monotonicity flag: Values --> OFF - high-order reconstruction,
                                   //                                  ON - limited linear reconstruction
  int rings;                       //!< Number of rings used to generate the reconstruction stencil
  
  short PseudoInverseFlag;	   //!< Flag to indicate whether the pseudo-inverse has been already computed.

  //! Storage for the pseudo-inverse of the LHS term in the CENO reconstruction.
  DenseMatrix **CENO_LHS;
  //! Storage for the geometric weights used in CENO reconstruction.
  ColumnVector **CENO_Geometric_Weights;

  // Associate this reconstruction to a certain mesh block
  GeometryType* Geom;            //!< Pointer to 2D mesh block

  //! Set Geometry_Ptr to point to the same geometry as the current Geom pointer
  void set_geometry_pointer(GeometryType* & Geometry_Ptr) const { Geometry_Ptr = Geom; }

  //! Return the required number of neighbour rings
  static int NumberOfRings(int number_of_Taylor_derivatives);

};

/******************************************************
 * Initialize the HighOrder2D class static variables  *
 *****************************************************/


/****************************************************
 * Implement the HighOrder2D class member functions *
 ***************************************************/
//! Default Constructor
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(void):TD(0), CENO_LHS(),
					   CENO_Geometric_Weights(), PseudoInverseFlag(OFF), Geom(NULL){
  // set rings
  SetRings();

  // initialize smoothness indicator
  InitializeSmoothnessIndicator();

  // initialize monotonicity flag
  InitializeMonotonicityFlag();

}

//! Main Constructor
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(int ReconstructionOrder, GeometryType & Block):
  TD(ReconstructionOrder), CENO_LHS(), CENO_Geometric_Weights(), PseudoInverseFlag(OFF){

  // set geometry pointer
  SetGeometryPointer(Block);

  // set rings
  SetRings();
  
  // initialize smoothness indicator
  InitializeSmoothnessIndicator();

  // initialize monotonicity flag
  InitializeMonotonicityFlag();

}

//! Copy constructor 
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(const HighOrder2D<SOLN_STATE> & rhs): Geom(NULL)
{
  TD = rhs.TaylorDeriv();
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (!rhs.GeomWeights().null()) ){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
    PseudoInverseFlag = ON;
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);
}

//! Assignment operator
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE> & HighOrder2D<SOLN_STATE>::operator=(const HighOrder2D<SOLN_STATE> & rhs){

  // Handle self-assignment:
  if (this == & rhs) return *this;

  TD = rhs.TaylorDeriv();
  SI = rhs.CellSmoothnessIndicator();
  LimitedCell = rhs.CellInadequateFit();
  rings = rhs.CellRings();

  if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (!rhs.GeomWeights().null()) ){
    CENO_LHS = rhs.LHS();
    CENO_Geometric_Weights = rhs.GeomWeights();
    PseudoInverseFlag = ON;
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);

  return *this;
}

// deallocate()
/*! Deallocate memory when high-order is not required.
 *  Automatic deallocation is already provided when the 
 *  variables are deleted.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::deallocate(void){
  // deallocate TD
  TaylorDeriv().free_memory();

  // deallocate LimitedCell
  CellInadequateFit().reserve(0);
}

// AssociateGeometry()
/*! Assign a specific geometry to the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::AssociateGeometry(GeometryType & Block){
  // Set geometry pointer
  SetGeometryPointer(Block);
}

// NumberOfRings()
/*! Return the required number of neighbor rings
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::NumberOfRings(int number_of_Taylor_derivatives){
  switch(number_of_Taylor_derivatives){
  case 1: 
    return 0;	// piecewise constant
  case 2:
    return 1; 	// piecewise linear
  case 3:
    return 2;	// piecewise quadratic
  case 4:
    return 2;	// piecewise cubic
  case 5:
    return 3;	// piecewise quartic
  case 6:
    return 3;	// piecewise quintic
  case 7:
    return 4;	// piecewise 6th-order
  case 8:
    return 4;	// piecewise 7th-order
  case 9:
    return 5;	// piecewise 8th-order
  default:
    return -1;
  }
}

// SetRings()
/*! Set the number of rings around the current cell
 * which will be used to form the supporting stencil
 * of the reconstruction. 
 * In 1D, each ring adds two neighbours to the stencil.
 * TD.size() --> gives the number of TaylorDerivatives,
 * which is a function of the order of reconstruction
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetRings(void){
  rings = NumberOfRings(TD.size());
}

//!< StencilSize()
/*!
 * Return the stencil size for CENO reconstruction.
 */
template<class SOLN_STATE> inline
int StencilSize(void) const {
  static int temp;

  // Calculate temp based on the number of rings
  temp = 1 + 2*CellRings();

  return temp*temp;
} 

// Nghost(ReconstructionOrder)
/*! Returns the minimum number of ghost cells for the solution domain
 *  required to carry out a reconstruction of order ReconstructionOrder
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::Nghost(int ReconstructionOrder){
  // Compute the size of the TaylorDerivative container
  int TD_size(ReconstructionOrder+1);

  // Return final number of ghost cells
  // --> double the number of rings: smoothness indicator requirement)
  // --> add one extra cell: boundary flux calculation requirement
  return 1 + 2* NumberOfRings(TD_size);
}

//! Reserve memory and set initial values for monotonicity flags.
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeMonotonicityFlag(void){
  LimitedCell.reserve(Soln_State::NumberOfVariables);
  for (int i = 0; i <= Soln_State::NumberOfVariables - 1; ++i){
    LimitedCell[i] = OFF; /* initialize the flags to OFF (smooth solution)*/
  }
}

//! Reserve memory and set initial values for smoothness indicators.
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeSmoothnessIndicator(void){
  SI.reserve(Soln_State::NumberOfVariables);
  for (int i = 0; i <= Soln_State::NumberOfVariables - 1; ++i){
    SI[i] = 0.0; /* initialize the smoothness indicator values */
  }
}

// InitializeVariable()
/*! 
 * Allocate memory and initialize the high-order variables
 * based on the ReconstructionOrder and the ReconstructionMethod.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::InitializeVariable(int ReconstructionOrder, int ReconstructionMethod){
  const int Reconstruction_Order(ReconstructionOrder);
  const int Reconstruction_Method(ReconstructionMethod);

  // Set specific variables to each reconstruction method
  switch(Reconstruction_Method){
  case RECONSTRUCTION_CENO:
    // Generate TaylorDerivatives container for the required reconstruction order
    TaylorDeriv().GenerateContainer(Reconstruction_Order);
    
    // Generate Geometric Coefficients container for the required reconstruction order
    //    CellGeomCoeff().GenerateContainer(Reconstruction_Order);
    
    // Set the number of rings
    SetRings();
    
    // Allocate memory for the pseudo-inverse and the vector of geometric weights
    if ( CENO_Execution_Mode::CENO_SPEED_EFFICIENT ){
      LHS().newsize(StencilSize() - 1, NumberOfTaylorDerivatives() - 1);
      GeomWeights().newsize(StencilSize() );
    }
    break;

  default:
    // Lower-order schemes
    deallocate();

    // Set the number of rings
    SetRings();
  }
}

// InitializeVariable()
/*! 
 * Initialize the high-order variables with the parameters provided by
 * an input parameters object.
 */
template<class SOLN_STATE>
template<typename InputParametersType> inline
void HighOrder2D<SOLN_STATE>::InitializeVariable(const InputParametersType & IP){
  InitializeVariable(IP.ReconstructionOrder(), IP.i_ReconstructionMethod);
}

//! Integrate over the cell
template<class SOLN_STATE> template<typename FO, class ReturnType> inline
ReturnType HighOrder2D<SOLN_STATE>::IntegrateOverTheCell(const FO FuncObj,
							 const int & digits,
							 ReturnType _dummy_param){
  return AdaptiveGaussianQuadrature(FuncObj,
				    CellCenter() - 0.5* CellDelta(),
				    CellCenter() + 0.5* CellDelta(),
				    _dummy_param,digits);
}

// ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::
ComputeUnlimitedSolutionReconstruction(Soln_Block_Type *SolnBlk,
				       const int iCell,
				       const int ReconstructionMethod,
				       HighOrder2D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void)){

  vector<int> i_index(StencilSize()); 
  string msg;

  switch(ReconstructionMethod){
  case RECONSTRUCTION_CENO:
    // Make Stencil
    MakeReconstructionStencil(CellRings(),iCell,i_index);

    // Solve reconstruction for the current cell
    kExact_Reconstruction(*this,SolnBlk,AccessToHighOrderVar,iCell,i_index,StencilSize(),NumberOfTaylorDerivatives());

    // Reset the CellInadequateFit flag & the limiter value in the Taylor derivatives container
    for (int i = 0; i <= Soln_State::NumberOfVariables - 1; ++i){
      LimitedCell[i] = OFF; /* reset the flags to OFF (smooth solution)*/
    }
    TaylorDeriv().ResetLimiter();
    
    break;
    
  default:
    throw runtime_error("HighOrder2D ERROR: Unknown specified reconstruction method");

  }
}

// ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeReconstructionPseudoInverse(Soln_Block_Type *SolnBlk,
								 const int iCell,
								 HighOrder2D<Soln_State> & 
								 (Soln_Block_Type::*AccessToHighOrderVar)(void)){

  if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
      CENO_Execution_Mode::CENO_SPEED_EFFICIENT && 
      PseudoInverseFlag == OFF){

    // == Check if the reconstruction polynomial is piecewise constant
    if (NumberOfTaylorDerivatives() == 1){
      // Set properly the PseudoInverseFlag
      PseudoInverseFlag = ON;
      return;
    }

    vector<int> i_index(StencilSize()); 

    // Make Stencil
    MakeReconstructionStencil(CellRings(),iCell,i_index);
    
    // Form the left-hand-side (LHS) term for the current cell
    kExact_Reconstruction_Compute_LHS(*this,SolnBlk,AccessToHighOrderVar,
				      iCell,i_index,StencilSize(),NumberOfTaylorDerivatives());
    
    // Compute the pseudo-inverse and override the LHS term
    LHS().pseudo_inverse_override();

    // Reset the CellInadequateFit flag & the limiter value in the Taylor derivatives container
    for (int i = 0; i <= Soln_State::NumberOfVariables - 1; ++i){
      LimitedCell[i] = OFF; /* reset the flags to OFF (smooth solution)*/
    }
    TaylorDeriv().ResetLimiter();

    // Set properly the PseudoInverseFlag
    PseudoInverseFlag = ON;
  }
}


/*! 
 * This reconstruction is carried out when the order or reconstruction
 * is required to be dropped since the high-order interpolant is detected
 * to be non-smooth. \n
 * The high-order interpolant is going to be overwritten by the low-order one.
 * \param [in] SolnBlk The solution domain
 * \param [in] iCell  The cell for which the reconstruction is done
 * \param [in] Limiter The limiter used during this linear least-squares reconstruction
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::ComputeLowOrderReconstruction(Soln_Block_Type *SolnBlk,
							    const int iCell,
							    const int Limiter){

  // Local variables
  int i, n, n2, n_pts, index[2];
  double u0Min, u0Max, uQuad[2], phi;
  double Dx, DxDx_ave;
  Soln_State DU, DUDx_ave, dWdx;
  int TD;

  /* Carry out the limited linear least-squares solution reconstruction. */

  n_pts = 2;
  index[0] = iCell-1;
  index[1] = iCell+1; 
    
  DUDx_ave = Soln_State(0);
  DxDx_ave = ZERO;
    
  for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
    Dx = SolnBlk[ index[n2] ].CellCenter() - SolnBlk[iCell].CellCenter();
    DU = SolnBlk[ index[n2] ].CellSolutionPrimVar() - SolnBlk[iCell].CellSolutionPrimVar();
    DUDx_ave += DU*Dx;
    DxDx_ave += Dx*Dx;
  } /* endfor */
    					    
  DUDx_ave = DUDx_ave/double(n_pts);
  DxDx_ave = DxDx_ave/double(n_pts);
	
  dWdx = DUDx_ave/DxDx_ave;
	
  for ( n = 1 ; n <= Soln_State::NumberOfVariables ; ++n ) {

    if (CellInadequateFit(n) == ON){ // drop the order only for the variables that are flagged as unfit

      /* Zero all the derivatives but the first two ones associated with this parameter. */
      for (TD = 2; TD<NumberOfTaylorDerivatives(); ++TD){
	TaylorDeriv(TD,n) = 0.0;
      }

      /* Copy the first order derivative in the derivatives container. */
      TaylorDeriv(0,n) = SolnBlk[iCell].CellSolutionPrimVar(n);
      TaylorDeriv(1,n) = dWdx[n];

      /* Compute the limiter value for this parameter */
      u0Min = SolnBlk[iCell].CellSolutionPrimVar(n);
      u0Max = u0Min;
      for ( n2 = 0 ; n2 <= n_pts-1 ; ++n2 ) {
	u0Min = min(u0Min, SolnBlk[ index[n2] ].CellSolutionPrimVar(n));
	u0Max = max(u0Max, SolnBlk[ index[n2] ].CellSolutionPrimVar(n));
      } /* endfor */

      uQuad[0] = SolnBlk[iCell].CellSolutionPrimVar(n) - HALF*dWdx[n]*SolnBlk[iCell].CellDelta();
      uQuad[1] = SolnBlk[iCell].CellSolutionPrimVar(n) + HALF*dWdx[n]*SolnBlk[iCell].CellDelta();

      switch(Limiter) {
      case LIMITER_BARTH_JESPERSEN :
	phi = Limiter_BarthJespersen(uQuad, SolnBlk[iCell].CellSolutionPrimVar(n), u0Min, u0Max, 2);
	break;
      case LIMITER_VENKATAKRISHNAN :
	phi = Limiter_Venkatakrishnan(uQuad, SolnBlk[iCell].CellSolutionPrimVar(n), u0Min, u0Max, 2);
	break;
      case LIMITER_VANLEER :
	phi = Limiter_VanLeer(uQuad, SolnBlk[iCell].CellSolutionPrimVar(n), u0Min, u0Max, 2);
	break;
      case LIMITER_VANALBADA :
	phi = Limiter_VanAlbada(uQuad, SolnBlk[iCell].CellSolutionPrimVar(n), u0Min, u0Max, 2);
	break;
      case LIMITER_ZERO :
	phi = ZERO;
	break;
      case LIMITER_ONE :
	phi = ONE;
	break;
      default:
	throw runtime_error("ComputeLowOrderReconstruction() ERROR: Unknown limiter type");
      } /* endswitch */

      /* Copy the limiter value to the derivatives container. */
      TaylorDeriv().Limiter(n) = phi;
    } // endif
  } /* endfor (n) */

}

// ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::ComputeSmoothnessIndicator(Soln_Block_Type *SolnBlk,
							 const int iCell,
							 HighOrder2D<Soln_State> &
							 (Soln_Block_Type::*AccessToHighOrderVar)(void)){

  static double SS_Regression, SS_Residual; // regression sum of squares, residual sum of squares
  static double MeanSolution, alpha;
  static int DOF; 		 /* degrees of freedom */
  static double AdjustmentCoeff; /* adjustment coefficient */
  static double Temp, DeltaTol;	
  static int parameter, cell, ComputeSI;
  int _StencilSize_(StencilSize());
  vector<int> i_index(_StencilSize_); 

  // Make Stencil
  MakeReconstructionStencil(CellRings(),iCell,i_index);

  /* Compute the CENO smoothness indicator for the current cell */

  /* Initialize the static variables */
  DOF = NumberOfTaylorDerivatives();
  AdjustmentCoeff = (_StencilSize_ - DOF)/(DOF - 1.0);

  for (parameter=1; parameter<=Soln_State::NumberOfVariables; ++parameter){

    // Assign the Mean Solution
    MeanSolution = SolnBlk[iCell].CellSolutionPrimVar(parameter);

    /* DeltaTolerance */
    DeltaTol = CENO_Tolerances::SquareToleranceAroundValue(MeanSolution);

    // Compute the regression and residual sums
    ComputeSI = OFF;		/* assume that the smoothness indicator is not computed but assigned */

    /* Initialize SS_Regression & SS_Residual with the values obtained for iCell */
    Temp = (SolnBlk[iCell].*AccessToHighOrderVar)().TaylorDeriv(0,parameter) - MeanSolution;
    Temp *= Temp;		/* compute Temp square */
    SS_Regression = Temp;
    
    /* Check if the Temp is greater than DeltaTol */
    if (Temp > DeltaTol){
      ComputeSI = ON; 	        /* Decide to compute the smoothness indicator */
    }
    SS_Residual = 0.0;		/* for iCell this term is 0.0 */
    
    /* compute S(quare)S(sum)_Regression and SS_Residual for the rest of the stencil */
    for(cell=1; cell<_StencilSize_; ++cell){

      Temp = (SolnBlk[i_index[cell]].*AccessToHighOrderVar)().TaylorDeriv(0,parameter) - MeanSolution;
      Temp *= Temp;		/* compute Temp square */

      if (CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS){
	if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT){
	  /* the weighting works only with the speed efficient CENO */
	  SS_Regression += (SolnBlk[iCell].*AccessToHighOrderVar)().GeomWeights(cell) * Temp;
	} else {
	  throw runtime_error("ComputeSmoothnessIndicator() 1D ERROR: CENO_CONSIDER_WEIGHTS only works with CENO_SPEED_EFFICIENT");
	}
      } else {
	SS_Regression += Temp;
      }
      
      /* Check if any of the Temps is greater than DeltaTol */
      if ((ComputeSI==OFF) && (Temp > DeltaTol)){
	ComputeSI = ON; 	/* Decide to compute the smoothness indicator */
      }
      
      Temp = ( (SolnBlk[i_index[cell]].*AccessToHighOrderVar)().TaylorDeriv(0,parameter) - 
	       (SolnBlk[iCell].*AccessToHighOrderVar)().SolutionAtCoordinates( (SolnBlk[i_index[cell]].*AccessToHighOrderVar)().CellCenter(),parameter) );

      if (CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS){
	if (CENO_Execution_Mode::CENO_SPEED_EFFICIENT){
	  /* the weighting works only with the speed efficient CENO */
	  SS_Residual += (SolnBlk[iCell].*AccessToHighOrderVar)().GeomWeights(cell) * Temp * Temp;
	}
      } else {
	SS_Residual += Temp * Temp;
      }
    }
    
    // Decide if the smoothness indicator is computed or not
    if (ComputeSI){ 
      alpha = 1.0 - SS_Residual/SS_Regression;
    } else {
      // Assign the perfect fit value to the smoothness indicator
      alpha = 1.0;
    }

    (SolnBlk[iCell].*AccessToHighOrderVar)().CellSmoothnessIndicator(parameter) = 
      (alpha/(max(CENO_Tolerances::epsilon,1.0 - alpha))) * AdjustmentCoeff;

  }//endfor -> parameter
}


//! Compute the solution at the left cell interface
template<class SOLN_STATE> inline
SOLN_STATE HighOrder2D<SOLN_STATE>::left_state(void){
  return SolutionAtCoordinates(CellCenter() - 0.5* CellDelta());
}

//! Compute the solution at the right cell interface
template<class SOLN_STATE> inline
SOLN_STATE HighOrder2D<SOLN_STATE>::right_state(void){
  return SolutionAtCoordinates(CellCenter() + 0.5* CellDelta());
}


/*! 
 * Compute the integral over the cell geometry of the error between the
 * reconstructed polynomial and the function provided as input. 
 *
 * \param [in] FuncObj  The function relative to which the error is evaluated
 * \param [in] parameter The parameter for which the reconstruction is evaluated
 */
template<class SOLN_STATE> 
template<typename Function_Object_Type>
double HighOrder2D<SOLN_STATE>::ComputeSolutionErrorL1(const Function_Object_Type FuncObj,
						       const unsigned parameter){
  
  // Set the type of the returned value
  double _dummy_param(0.0);

  // set the pointer to the member function that is used to compute the solution of the polynomial reconstruction
  MemberFunction_TwoArguments_OneParameter_Type ReconstructedSolution = &ClassType::SolutionAtCoordinates;

  // Call the integration function
  return IntegrateOverTheCell(error_function(FuncObj,
					     wrapped_member_function_one_parameter(this,
										   ReconstructedSolution,
										   parameter,
										   _dummy_param),
					     _dummy_param),
			      10,_dummy_param);
}

template<class SOLN_STATE> 
double HighOrder2D<SOLN_STATE>::ComputeSolutionErrorL1(HighOrder2D<Soln_State> & Obj,
						       const unsigned parameter){

  // Set the type of the returned value
  double _dummy_param(0.0);

  // set the pointer to the member function that is used to compute the solution of the polynomial reconstruction
  MemberFunction_TwoArguments_OneParameter_Type ReconstructedSolution = &ClassType::SolutionAtCoordinates;

  // Call the computation of the error routine with a function given by a wrapper based on the argument
  return ComputeSolutionErrorL1(wrapped_member_function_one_parameter(&Obj,
								      ReconstructedSolution,
								      parameter,
								      _dummy_param),
				parameter);
}

/*! 
 * Compute the integral over the cell geometry of the squared error between the
 * reconstructed polynomial and the function provided as input. 
 *
 * \param [in] FuncObj  The function relative to which the error is evaluated
 * \param [in] parameter The parameter for which the reconstruction is evaluated
 */
template<class SOLN_STATE> 
template<typename Function_Object_Type>
double HighOrder2D<SOLN_STATE>::ComputeSolutionErrorL2(const Function_Object_Type FuncObj, const unsigned parameter){

  // Set the type of the returned value
  double _dummy_param(0.0);

  // set the pointer to the member function that is used to compute the solution of the polynomial reconstruction
  MemberFunction_TwoArguments_OneParameter_Type ReconstructedSolution = &ClassType::SolutionAtCoordinates;

  // Call the integration function
  return IntegrateOverTheCell(square_error_function(FuncObj,
						    wrapped_member_function_one_parameter(this,
											  ReconstructedSolution,
											  parameter,
											  _dummy_param),
						    _dummy_param),
			      10,_dummy_param);
}


// Friend functions
//! operator== 
template<class SOLN_STATE> inline
bool operator== (const HighOrder2D<SOLN_STATE> & left, const HighOrder2D<SOLN_STATE> & right){
   
  if ( CENO_Execution_Mode::CENO_SPEED_EFFICIENT && (!(left.GeomWeights().null() && right.GeomWeights().null())) ){
    return ( (left.TaylorDeriv() == right.TaylorDeriv()) && 
	     (left.CellSmoothnessIndicator() == right.CellSmoothnessIndicator()) &&
	     (left.CellInadequateFit() == right.CellInadequateFit()) &&
	     (left.CellRings() == right.CellRings()) && 
	     (left.CellGeometry() == right.CellGeometry()) && 
	     (left.LHS() == right.LHS()) &&
	     (left.GeomWeights() == right.GeomWeights()) 
	     );
  } else {
    // There is no memory allocated for pseudo-inverse data in both objects
    return ( (left.TaylorDeriv() == right.TaylorDeriv()) && 
	     (left.CellSmoothnessIndicator() == right.CellSmoothnessIndicator()) &&
	     (left.CellInadequateFit() == right.CellInadequateFit()) &&
	     (left.CellRings() == right.CellRings()) && 
	     (left.CellGeometry() == right.CellGeometry()) );
  }
}

//! operator!= 
template<class SOLN_STATE> inline
bool operator!= (const HighOrder2D<SOLN_STATE> & left, const HighOrder2D<SOLN_STATE> & right){
  return !(left == right);
}

//! operator<<
template<class SOLN_STATE> inline
ostream & operator<< (ostream & os, const HighOrder2D<SOLN_STATE> & Obj){

  os.setf(ios::skipws,ios::scientific);
  os << Obj.TaylorDeriv();
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
istream & operator>> (istream & os, HighOrder2D<SOLN_STATE> & Obj){

  os.setf(ios::skipws);
  os >> Obj.TaylorDeriv()
     >> Obj.CellSmoothnessIndicator();

  // Set rings
  Obj.SetRings();
  os.unsetf(ios::skipws);
  return os;
}


// Specializations
template<> inline
void HighOrder2D<double>::InitializeMonotonicityFlag(void){
  LimitedCell.reserve(1);
  LimitedCell[0] = OFF;
}

template<> inline
double & HighOrder2D<double>::TaylorDeriv(const int & p1, const int & Variable) {
 return TD(p1);
}

template<> inline
const double & HighOrder2D<double>::TaylorDeriv(const int & p1, const int & Variable) const {
 return TD(p1);
}

template<> inline
const double & HighOrder2D<double>::CellSmoothnessIndicator(const int & VarPosition) const {
  return SI;
}

template<> inline
double & HighOrder2D<double>::CellSmoothnessIndicator(const int & VarPosition) {
  return SI;
}

template <> inline
double HighOrder2D<double>::SolutionAtCoordinates(const double & X_Coord, const unsigned parameter){
  return TD.ComputeSolutionFor(X_Coord - CellCenter());
}



// HighOrderSolutionReconstructionOverDomain()
/*! 
 * Compute the high-order reconstruction for each computational cell 
 * of the SolnBlk using the 'IP.i_ReconstructionMethod' algorithm.
 *
 * \param IP input parameter object. Provides the reconstruction method
 * \param AccessToHighOrderVar member function of Soln_Block_Type 
 * that returns the high-order variable which is used in the
 * reconstruction process.
 */
template<class Soln_Block_Type, class InputParametersType>
void HighOrderSolutionReconstructionOverDomain(Soln_Block_Type *SolnBlk,
					       const InputParametersType & IP,
					       typename Soln_Block_Type::HighOrderType & 
					       (Soln_Block_Type::*AccessToHighOrderVar)(void)) {

  typedef typename Soln_Block_Type::HighOrderType HighOrderType;

  int ICl(SolnBlk[0].ICl), ICu( SolnBlk[0].ICu);
  int i, parameter;
  bool InadequateFitFlag;

  switch(IP.i_ReconstructionMethod){
    /* C(entral)ENO -> central stencil with post-analysis of the reconstruction */
  case RECONSTRUCTION_CENO:
    // require a minimum number of ghost cells equal to what is necessary for the current high-order reconstruction
    require(SolnBlk[0].Nghost >= HighOrderType::Nghost((SolnBlk[0].*AccessToHighOrderVar)().CellRecOrder()),
	    "ReconstructSolutionOverDomain() ERROR: Not enough ghost cells to perform the current reconstruction");

    //Step 1: Compute the k-exact reconstruction
    for (i = ICl - ((SolnBlk[0].*AccessToHighOrderVar)().CellRings() + 1);
	 i<= ICu + ((SolnBlk[0].*AccessToHighOrderVar)().CellRings() + 1);
	 ++i) {

      // Compute PseudoInverse if required
      (SolnBlk[i].*AccessToHighOrderVar)().ComputeReconstructionPseudoInverse(SolnBlk,i);

      // Compute Unlimited High-Order Reconstruction
      (SolnBlk[i].*AccessToHighOrderVar)().ComputeUnlimitedSolutionReconstruction(SolnBlk,i,RECONSTRUCTION_CENO,
										  AccessToHighOrderVar);
    }
    
    // Step 2 and 3: Check smoothness
    for (i=ICl-1; i<=ICu+1; ++i){
      
      //Step 2: Compute the Smoothness Indicator for the cells used to compute the Riemann problem.
      (SolnBlk[i].*AccessToHighOrderVar)().ComputeSmoothnessIndicator(SolnBlk,i,AccessToHighOrderVar);
      
      //Step 3: Do a post-reconstruction analysis
      /* Check the smoothness condition */
      for(parameter=1; parameter<=Soln_Block_Type::HighOrderType::Soln_State::NumberOfVariables; ++parameter){
	if( (SolnBlk[i].*AccessToHighOrderVar)().CellSmoothnessIndicator(parameter) < CENO_Tolerances::Fit_Tolerance ){

	  /* Flag the 'i' cell with non-smooth reconstruction */
	  (SolnBlk[i].*AccessToHighOrderVar)().CellInadequateFit(parameter) = ON;

	  if (CENO_Execution_Mode::CENO_PADDING){
	    /* Flag all the cell surrounding the 'i' cell with bad reconstruction if CENO_Padding is ON */
	    (SolnBlk[i-1].*AccessToHighOrderVar)().CellInadequateFit(parameter) = ON;
	    (SolnBlk[i+1].*AccessToHighOrderVar)().CellInadequateFit(parameter) = ON;
	  }
	}//endif
      }//endfor(parameter)
      
    } //endfor(i)
    
    //Step 4: Switch the high-order reconstruction to a monotone piecewise one for 
    //        those cells that are detected as unfit.
    for (i=ICl-1; i<=ICu+1; ++i){
      
      // Reset flag
      InadequateFitFlag = false;
      
      // analyse the 'CellInadequateFit' flags and set 'InadequateFitFlag'
      for(parameter=1; parameter<=Soln_Block_Type::HighOrderType::Soln_State::NumberOfVariables; ++parameter){
	if (InadequateFitFlag == true){	// break the loop if the flag is already 'true'
	  break;
	} else if ( (SolnBlk[i].*AccessToHighOrderVar)().CellInadequateFit(parameter) == ON ){
	  InadequateFitFlag = true;
	}
      }//endfor(parameter)
      
      if (InadequateFitFlag == true && CENO_Execution_Mode::CENO_DROP_ORDER){
	(SolnBlk[i].*AccessToHighOrderVar)().ComputeLowOrderReconstruction(SolnBlk,i,IP.i_Limiter);
      }

    }//endfor (i) 

    break;
    
  default:
    throw runtime_error("ReconstructSolutionOverDomain ERROR: Unknown reconstruction method!");
  } /* endswitch */
  
}

#endif
