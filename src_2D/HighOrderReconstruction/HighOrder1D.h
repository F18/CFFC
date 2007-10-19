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
#include "ReconstructionSolvers1D.h"

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
 * \nosubgrouping
 *********************************************************/
template<class SOLN_STATE>
class HighOrder1D{
  
public: 

  //@{ @name Defined public types:
  typedef SOLN_STATE Soln_State;
  typedef Cell1D_Uniform GeometryType;
  typedef TaylorDerivativesContainer<OneD,Soln_State> DerivativesContainer;
  typedef TaylorDerivativesContainer<OneD,double> GeometricMoments;
  typedef typename DerivativesContainer::Derivative  Derivative;
  //@}
  
  //@{ @name Constructors:
  HighOrder1D(void);		//!< Simple constructor 
  HighOrder1D(int ReconstructionOrder, GeometryType & Cell); //!< Advanced constructor
  HighOrder1D( const HighOrder1D & rhs); //!< Copy constructor
  //@}


  //@{ @name Destructors:
  ~HighOrder1D(void){ Geom = NULL;}
  void deallocate(void);
  //@} 

  HighOrder1D<Soln_State> & operator=(const HighOrder1D<Soln_State> & rhs); //!< Assignment operator

  //@{ @name Field access functions:
  //@{ @name Taylor Derivatives:
  const DerivativesContainer & CellDeriv(void) const {return TD;}
  DerivativesContainer & CellDeriv(void) {return TD;}
  const Soln_State & CellDeriv(const int & p1) const {return TD(p1);}
  Soln_State & CellDeriv(const int & p1) {return TD(p1);}
  const double & CellDeriv(const int & p1, const int & Variable) const { return TD(p1)[Variable];}
  double & CellDeriv(const int & p1, const int & Variable) { return TD(p1)[Variable];}
  Derivative & CellDeriv_InPosition(const int & position) {return TD(position,true,true,true);}
  const Derivative & CellDeriv_InPosition(const int & position) const {return TD(position,true,true,true);}
  //@} 

  const int NumberOfTaylorDerivatives(void) const {return TD.size();}

  const GeometricMoments & CellGeomCoeff(void) const { return GeomCoeff; }
  GeometricMoments & CellGeomCoeff(void) {return GeomCoeff;}
  const double & CellGeomCoeff(const int & p1) {return GeomCoeff(p1);}
  const double & CellGeomCoeff_InPosition(const int & position) {return GeomCoeff(position,true,true,true).D();}

  const int & CellRings(void) const {return rings;}
  const int CellRecOrder(void) const {return TD.RecOrder();}

  MemoryStorageENO_1D & ENO_MemoryPool(void){ return ENO_Mem; }
  const MemoryStorageENO_1D & ENO_MemoryPool(void) const { return ENO_Mem; }

  //@{ @name Monotonicity info for high-order
  const vector<short int> & CellInadequateFit(void) const { return LimitedCell;}
  vector<short int> & CellInadequateFit(void){ return LimitedCell;}
  const short int & CellInadequateFit( const int VarPosition) const { return LimitedCell[VarPosition-1];}
  short int & CellInadequateFit( const int VarPosition){ return LimitedCell[VarPosition-1];}
  //@}

  //@{ @name Smoothness indicator
  const Soln_State & CellSmoothnessIndicator(void) const { return SI;}
  Soln_State & CellSmoothnessIndicator(void){ return SI;}
  const double & CellSmoothnessIndicator(const int & VarPosition)const{ return SI[VarPosition];}
  double & CellSmoothnessIndicator( const int & VarPosition){ return SI[VarPosition];}
  //@}

  //@{ @name Pseudo-inverse of the LHS term in the CENO reconstruction
  DenseMatrix & LHS(void) {return CENO_LHS;}
  const DenseMatrix & LHS(void) const {return CENO_LHS;}
  double & LHS(const int & IndexI, const int & IndexJ) {return CENO_LHS(IndexI,IndexJ);}
  const double & LHS(const int & IndexI, const int & IndexJ) const {return CENO_LHS(IndexI,IndexJ);}
  //@}

  //@{ @name Geometric weights assigned to the cells that are part of the stencil
  ColumnVector & GeomWeights(void){return CENO_Geometric_Weights;}
  const ColumnVector & GeomWeights(void) const {return CENO_Geometric_Weights;}
  double & GeomWeights(const int & Index){return CENO_Geometric_Weights(Index);}
  const double & GeomWeights(const int & Index) const {return CENO_Geometric_Weights(Index);}
  //@}
  //@}

  //@{ @name Cell geometry
  GeometryType* Geometry(void){return Geom;}
  const GeometryType* Geometry(void) const {return Geom;}
  GeometryType & CellGeometry(void){return *Geom;}
  const GeometryType & CellGeometry(void) const {return *Geom;}
  const double & CellCenter(void) const {return Geom->x;}
  const double & CellDelta (void) {return Geom->dx;}
  void SetGeometryPointer(GeometryType & Cell){ Geom = &Cell;}
  void AssociateGeometry(GeometryType & Cell);
  //@}

  /* Operating functions */
  static int Nghost(int ReconstructionOrder); //!< return the required number of ghost cells
  void SetRings(void);
  int StencilSize(void) const {return 1 + 2*CellRings();} //!< return the stencil size for CENO reconstruction
  void ComputeGeometricCoefficients(void);
  void InitializeMonotonicityFlag(void);

  void InitializeVariable(int ReconstructionOrder, int ReconstructionMethod);
  template<typename InputParametersType>
  void InitializeVariable(const InputParametersType & IP);

  double SolutionAtCoordinates(const double & X_Coord, const unsigned parameter){
    return TD.ComputeSolutionFor(X_Coord - CellCenter())[parameter];
  }
  Soln_State SolutionAtCoordinates(const double & X_Coord){
    return TD.ComputeSolutionFor(X_Coord - CellCenter());
  }
  /* Computation of right and left states (left and right are relative to the position of the centroid) */
  Soln_State left_state(void);
  Soln_State right_state(void);

  /*!
   * Integrate over the domain of the geometry associated with this high-order solution
   */
  template<typename FO, class ReturnType>
  ReturnType IntegrateOverTheCell(const FO FuncObj, const int & digits, ReturnType _dummy_param);
  
  bool IsPseudoInversePreComputed(void) const { return (PseudoInverseFlag==ON) ? true:false; }
  void MustUpdatePseudoInverse(void) {PseudoInverseFlag = OFF;}

  //@{ @name Reconstructions:
  /*! 
   * Compute the unlimited high-order solution reconstruction.
   */
  template<class Soln_Block_Type>
  void ComputeUnlimitedSolutionReconstruction(Soln_Block_Type *SolnBlk,
					      const int iCell,
					      const int ReconstructionMethod,
					      HighOrder1D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
					      &Soln_Block_Type::CellHighOrder);

  /*! 
   * Compute the pseudo-inverse corresponding to the unlimited high-order solution reconstruction.
   */
  template<class Soln_Block_Type>
  void ComputeReconstructionPseudoInverse(Soln_Block_Type *SolnBlk,
					  const int iCell,
					  HighOrder1D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
					  &Soln_Block_Type::CellHighOrder);

  /*! 
   * Compute the second (low-order) reconstruction in the CENO algorithm.
   */
  template<class Soln_Block_Type>
  void ComputeLowOrderReconstruction(Soln_Block_Type *SolnBlk, const int iCell, const int Limiter);

  //@}

  //@{ @name CENO Analysis:
  template<class Soln_Block_Type>
  void ComputeSmoothnessIndicator(Soln_Block_Type *SolnBlk,
				  const int iCell,
				  HighOrder1D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void) = 
				  &Soln_Block_Type::CellHighOrder);
  //@}

  /* Friend functions */
  friend bool operator== <Soln_State> (const HighOrder1D<Soln_State> & left,
				       const HighOrder1D<Soln_State> & right);

  friend bool operator!= <Soln_State> (const HighOrder1D<Soln_State> & left,
				       const HighOrder1D<Soln_State> & right);

  friend ostream & operator<< <Soln_State> (ostream & os, const HighOrder1D<Soln_State> & Obj);
  friend istream & operator>> <Soln_State> (istream & os, HighOrder1D<Soln_State> & Obj);


protected:
  
private:
  DerivativesContainer TD;  	//!< High-order TaylorDerivatives
  GeometricMoments GeomCoeff;   //!< The integrals of the geometric moments with respect to the centroid
  Soln_State SI;                //!< The values of the smoothness indicator calculated for each variable
  vector<short int> LimitedCell; //!< Monotonicity flag: Values --> OFF - high-order reconstruction
                                 //                               ON - limited linear reconstruction
  int rings;                    //!< Number of rings used to generate the reconstruction stencil
  

  short PseudoInverseFlag;	//!< Flag to indicate whether the pseudo-inverse has been set or not

  /* Create storage for the pseudo-inverse of the LHS term in the CENO reconstruction */
  DenseMatrix CENO_LHS;
  ColumnVector CENO_Geometric_Weights;

  // Associate this reconstruction to a certain cell
  GeometryType* Geom;    //!< Pointer to cell geometry

  /* Create memory pool for the ENO reconstruction */
  static MemoryStorageENO_1D ENO_Mem;

  // Set Geometry_Ptr to point to the same geometry as the current Geom pointer
  void set_geometry_pointer(GeometryType* & Geometry_Ptr) const { Geometry_Ptr = Geom; }

  // Return the required number of neighbor rings
  static int NumberOfRings(int number_of_Taylor_derivatives);

};

/******************************************************
 * Initialize the HighOrder1D class static variables  *
 *****************************************************/
template<class SOLN_STATE> MemoryStorageENO_1D HighOrder1D<SOLN_STATE>::ENO_Mem = MemoryStorageENO_1D(0);

/****************************************************
 * Implement the HighOrder1D class member functions *
 ***************************************************/
//! Default Constructor
template<class SOLN_STATE> inline
HighOrder1D<SOLN_STATE>::HighOrder1D(void):TD(0), GeomCoeff(0), SI(0),
					   CENO_LHS(), CENO_Geometric_Weights(), PseudoInverseFlag(OFF), Geom(NULL) {
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
  SI(0), CENO_LHS(), CENO_Geometric_Weights(), PseudoInverseFlag(OFF){

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
    PseudoInverseFlag = ON;
  }

  // point to the same geometry as rhs.Geom
  rhs.set_geometry_pointer(Geom);

  return *this;
}

//! deallocate()
/*! Deallocate memory when high-order is not required.
 *  Automatic deallocation is already provided when the 
 *  variables are deleted.
 */
template<class SOLN_STATE> inline
void HighOrder1D<SOLN_STATE>::deallocate(void){
  // deallocate TD
  CellDeriv().free_memory();

  // deallocate geometric coefficients
  CellGeomCoeff().free_memory();

  // deallocate LimitedCell
  CellInadequateFit().reserve(0);
}

//! AssociateGeometry()
/*! Assign a specific geometry to the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder1D<SOLN_STATE>::AssociateGeometry(GeometryType & Cell){
  // Set geometry pointer
  SetGeometryPointer(Cell);

  // Recalculate the geometric coefficients
  ComputeGeometricCoefficients(); 
}

//! NumberOfRings()
/*! Return the required number of neighbor rings
 */
template<class SOLN_STATE> inline
int HighOrder1D<SOLN_STATE>::NumberOfRings(int number_of_Taylor_derivatives){
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
  rings = NumberOfRings(TD.size());
}

//! Nghost(ReconstructionOrder)
/*! Returns the minimum number of ghost cells for the solution domain
 *  required to carry out a reconstruction of order ReconstructionOrder
 */
template<class SOLN_STATE> inline
int HighOrder1D<SOLN_STATE>::Nghost(int ReconstructionOrder){
  // Compute the size of the TaylorDerivative container
  int TD_size(ReconstructionOrder+1);

  // Return final number of ghost cells
  // --> double the number of rings: smoothness indicator requirement)
  // --> add one extra cell: boundary flux calculation requirement
  return 1 + 2* NumberOfRings(TD_size);
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

  if (GeomCoeff.RecOrder() > -1){ // check for allocated memory
    /* Create the polynomial function for the cell */
    GeneralizedPolynomialFunctionOfOneVariable Polynom(0,CellCenter());
    GeomCoeff(0) = 1.0;
    for (int p1=1; p1<=GeomCoeff.RecOrder(); ++p1){
      Polynom.ChangePowerTo(p1);
      GeomCoeff(p1) = IntegrateOverTheCell(Polynom,14,GeomCoeff(p1))/CellDelta();
    }
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

//! InitializeVariable()
/*! 
 * Allocate memory and initialize the high-order variables
 * based on the ReconstructionOrder and the ReconstructionMethod.
 * 
 */
template<class SOLN_STATE>
void HighOrder1D<SOLN_STATE>::InitializeVariable(int ReconstructionOrder, int ReconstructionMethod){
  const int Reconstruction_Order(ReconstructionOrder);
  const int Reconstruction_Method(ReconstructionMethod);

  // Set specific variables to each reconstruction method
  switch(Reconstruction_Method){
  case RECONSTRUCTION_CENO:
    // Generate TaylorDerivatives container for the required reconstruction order
    CellDeriv().GenerateContainer(Reconstruction_Order);
    
    // Generate Geometric Coefficients container for the required reconstruction order
    CellGeomCoeff().GenerateContainer(Reconstruction_Order);
    
    // Set the number of rings
    SetRings();
    
    // Allocate memory for the pseudo-inverse and the vector of geometric weights
    if ( CENO_Execution_Mode::CENO_SPEED_EFFICIENT ){
      LHS().newsize(StencilSize() - 1, NumberOfTaylorDerivatives() - 1);
      GeomWeights().newsize(StencilSize() );
    }
    break;

  case RECONSTRUCTION_ENO:
    // Generate TaylorDerivatives container for the required reconstruction order
    CellDeriv().GenerateContainer(Reconstruction_Order);

    // Generate Geometric Coefficients container for the required reconstruction order
    CellGeomCoeff().GenerateContainer(Reconstruction_Order);

    // deallocate LimitedCell
    CellInadequateFit().reserve(0);

    // Set the number of rings
    SetRings();

    // Allocate enough memory in the memory pool
    if (ENO_Mem.NumOfUnknowns != CellDeriv().size()){
      ENO_Mem.newsize(CellDeriv().size());
    }
    break;

  case RECONSTRUCTION_ENO_CHARACTERISTIC:
    // Generate TaylorDerivatives container for the required reconstruction order
    CellDeriv().GenerateContainer(Reconstruction_Order);

    // Generate Geometric Coefficients container for the required reconstruction order
    CellGeomCoeff().GenerateContainer(Reconstruction_Order);

    // deallocate LimitedCell
    CellInadequateFit().reserve(0);

    // Set the number of rings
    SetRings();

    // Allocate enough memory in the memory pool
    if (ENO_Mem.NumOfUnknowns != CellDeriv().size()){
      ENO_Mem.newsize(CellDeriv().size());
    }
    break;

  default:
    // Lower-order schemes
    deallocate();

    // Set the number of rings
    SetRings();
  }
}

//! InitializeVariable()
/*! 
 * Initialize the high-order variables with the parameters provided by
 * an input parameters object.
 */
template<class SOLN_STATE>
template<typename InputParametersType> inline
void HighOrder1D<SOLN_STATE>::InitializeVariable(const InputParametersType & IP){
  InitializeVariable(IP.ReconstructionOrder(), IP.i_ReconstructionMethod);
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

//! ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder1D<SOLN_STATE>::
ComputeUnlimitedSolutionReconstruction(Soln_Block_Type *SolnBlk,
				       const int iCell,
				       const int ReconstructionMethod,
				       HighOrder1D<Soln_State> & (Soln_Block_Type::*AccessToHighOrderVar)(void)){

  vector<int> i_index(StencilSize()); 

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
    CellDeriv().ResetLimiter();
    
    break;
    
  case RECONSTRUCTION_ENO:
    ENO_Reconstruction(*this,SolnBlk,AccessToHighOrderVar,iCell,ENO_MemoryPool());
    break;

  case RECONSTRUCTION_ENO_CHARACTERISTIC:
    ENO_Characteristics_Reconstruction(*this,SolnBlk,AccessToHighOrderVar,iCell,ENO_MemoryPool());
    break;
    
  default:
    throw runtime_error("HighOrder1D ERROR: Unknown specified reconstruction method");

  }
}

//! ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder1D<SOLN_STATE>::ComputeReconstructionPseudoInverse(Soln_Block_Type *SolnBlk,
								 const int iCell,
								 HighOrder1D<Soln_State> & 
								 (Soln_Block_Type::*AccessToHighOrderVar)(void)){

  if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
      CENO_Execution_Mode::CENO_SPEED_EFFICIENT && 
      PseudoInverseFlag == OFF){

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
    CellDeriv().ResetLimiter();

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
void HighOrder1D<SOLN_STATE>::ComputeLowOrderReconstruction(Soln_Block_Type *SolnBlk,
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
	CellDeriv(TD,n) = 0.0;
      }

      /* Copy the first order derivative in the derivatives container. */
      CellDeriv(0,n) = SolnBlk[iCell].CellSolutionPrimVar(n);
      CellDeriv(1,n) = dWdx[n];

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
      CellDeriv().Limiter(n) = phi;
    } // endif
  } /* endfor (n) */

}

//! ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder1D<SOLN_STATE>::ComputeSmoothnessIndicator(Soln_Block_Type *SolnBlk,
							 const int iCell,
							 HighOrder1D<Soln_State> &
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
    DeltaTol = CENO_EpsilonTol::SquareToleranceAroundValue(MeanSolution);

    // Compute the regression and residual sums
    ComputeSI = OFF;		/* assume that the smoothness indicator is not computed but assigned */

    /* Initialize SS_Regression & SS_Residual with the values obtained for iCell */
    Temp = (SolnBlk[iCell].*AccessToHighOrderVar)().CellDeriv(0,parameter) - MeanSolution;
    Temp *= Temp;		/* compute Temp square */
    SS_Regression = Temp;
    
    /* Check if the Temp is greater than DeltaTol */
    if (Temp > DeltaTol){
      ComputeSI = ON; 	        /* Decide to compute the smoothness indicator */
    }
    SS_Residual = 0.0;		/* for iCell this term is 0.0 */
    
    /* compute S(quare)S(sum)_Regression and SS_Residual for the rest of the stencil */
    for(cell=1; cell<_StencilSize_; ++cell){

      Temp = (SolnBlk[i_index[cell]].*AccessToHighOrderVar)().CellDeriv(0,parameter) - MeanSolution;
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
      
      Temp = ( (SolnBlk[i_index[cell]].*AccessToHighOrderVar)().CellDeriv(0,parameter) - 
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
      (alpha/(max(CENO_EpsilonTol::epsilon,1.0 - alpha))) * AdjustmentCoeff;

  }//endfor -> parameter
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



//! ReconstructSolutionOverDomain()
/*! 
 * Compute the high-order reconstruction for each computational cell 
 * of the SolnBlk using the 'IP.i_ReconstructionMethod' algorithm.
 * 'AccessToHighOrderVar' is the member function of Soln_Block_Type 
 * that returns the high-order variable which is used in the
 * reconstruction process.
 */
template<class Soln_Block_Type, class InputParametersType>
void ReconstructSolutionOverDomain(Soln_Block_Type *SolnBlk,
				   const InputParametersType & IP,
				   HighOrder1D<typename Soln_Block_Type::HighOrderType::Soln_State> & 
				   (Soln_Block_Type::*AccessToHighOrderVar)(void)) {

  typedef HighOrder1D<typename Soln_Block_Type::HighOrderType::Soln_State> HighOrderType;

  int ICl(SolnBlk[0].ICl), ICu( SolnBlk[0].ICu);
  int i, parameter;
  bool InadequateFitFlag;

  switch(IP.i_ReconstructionMethod){
    /* C(entral)ENO -> central stencil with post-analysis of the reconstruction */
  case RECONSTRUCTION_CENO:
    // require a minimum number of ghost cells equal to what is necessary for the current high-order reconstruction
    require(SolnBlk[0].Nghost >= HighOrderType::Nghost((SolnBlk[i].*AccessToHighOrderVar)().CellRecOrder()),
	    "ReconstructSolutionOverDomain() ERROR: Not enough ghost cells to perform the current reconstruction");

    //Step 1: Compute the k-exact reconstruction
    for (i = ICl - ((SolnBlk[i].*AccessToHighOrderVar)().CellRings() + 1);
	 i<= ICu + ((SolnBlk[i].*AccessToHighOrderVar)().CellRings() + 1);
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
	if( (SolnBlk[i].*AccessToHighOrderVar)().CellSmoothnessIndicator(parameter) < IP.FitTolerance() ){

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
    
    /* Essentially NonOscillatory (ENO) */
  case RECONSTRUCTION_ENO:
    for (i = ICl-1 ;i<= ICu+1 ;++i) {
      /* Solve reconstruction for the current cell */
      (SolnBlk[i].*AccessToHighOrderVar)().ComputeUnlimitedSolutionReconstruction(SolnBlk,i,RECONSTRUCTION_ENO,
										  AccessToHighOrderVar);
    }
    break;
    
    /* Essentially NonOscillatory (ENO) in characteristic variables */
  case RECONSTRUCTION_ENO_CHARACTERISTIC:
    for (i = ICl-1 ;i<= ICu+1 ;++i) {
      /* Solve reconstruction for the current cell */
      (SolnBlk[i].*AccessToHighOrderVar)().ComputeUnlimitedSolutionReconstruction(SolnBlk,i,RECONSTRUCTION_ENO_CHARACTERISTIC,
										  AccessToHighOrderVar);
    }
    break;
    
  default:
    throw runtime_error("ReconstructSolutionOverDomain ERROR: Unknown reconstruction method!");
  } /* endswitch */
  
}

#endif
