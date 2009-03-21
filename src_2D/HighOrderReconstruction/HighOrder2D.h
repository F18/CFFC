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
#include "../Grid/HO_Grid2DQuad.h"
#include "ReconstructionHelpers.h"
#include "Cauchy_BoundaryConditions.h"
#include "HighOrder2D_BlockBoundary.h"

/*********************************
 * Declare the HighOrder2D class *
 ********************************/
template <class SOLN_STATE> 
class HighOrder2D;

/************************************************
 *     Friend Functions : HighOrder2D           *
 ************************************************/
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
  typedef std::vector<int> FlagType;
  //! Type for smoothness indicator associated with each reconstructed variable.
  typedef std::vector<double> DoubleArrayType;
  //! Type for array of Vector2D
  typedef vector<Vector2D> Vector2DArray;
  //! Type for high-order Cauchy boundary conditions
  typedef Cauchy_BCs<Soln_State> BC_Type;
  //! Type for the array of high-order solution boundary conditions
  typedef vector<BC_Type *> BC_Type_Array;
  //@}
  
  //! @name Constructors:
  //@{
  HighOrder2D(void);		//!< Simple constructor 
  HighOrder2D(int ReconstructionOrder, GeometryType & Block,
	      const bool &_pseudo_inverse_allocation_ = false); //!< Advanced constructor
  HighOrder2D( const HighOrder2D & rhs); //!< Copy constructor
  void allocate(const int &NC_IDir, const int &NC_JDir,
		const int &Nghost,
		const bool &_pseudo_inverse_allocation_,
		int ReconstructionOrder = -1);
  //@}
  

  //! @name Destructors:
  //@{
  ~HighOrder2D(void){ deallocate(); }
  void deallocate_CellMemory(void);
  void deallocate(void);
  void deallocate_ReconstructionTypeMap(void);
  //@} 

  HighOrder2D<Soln_State> & operator=(const HighOrder2D<Soln_State> & rhs); //!< Assignment operator

  //! @name Taylor derivatives:
  //@{
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
  //! Get the solution state of the Taylor derivative of cell (ii,jj) which is stored in the 'position' place.
  Soln_State & CellTaylorDerivState(const int & ii, const int & jj,
				    const int & position) {return TD[ii][jj](position).D();}
  //! Get the solution state of the Taylor derivative of cell (ii,jj) which is stored in the 'position' place.
  const Soln_State & CellTaylorDerivState(const int & ii, const int & jj,
					  const int & position) const {return TD[ii][jj](position).D();}
  //! Get the number of Taylor derivatives in each cell container
  const int NumberOfTaylorDerivatives(void) const {return TD[ICl][JCl].size();}
  //@} 

  //! @name Geometric moments
  //@{
  //! Get the geometric coefficients of cell (ii,jj).
  const GeometricMoments & CellGeomCoeff(const int & ii, const int & jj) const { return Geom->CellGeomCoeff(ii,jj); }
  //! Get the geometric coefficients of cell (ii,jj) which is stored in the 'position' place.
  const GeomMoment & CellGeomCoeff(const int & ii, const int & jj,
				   const int & position) {return Geom->CellGeomCoeff(ii,jj,position);}
  //! Get the value of cell (ii,jj) geometric coefficient with x-power 'p1' and y-power 'p2'.
  const double & CellGeomCoeffValue(const int & ii, const int & jj,
				    const int & p1, const int & p2) {return Geom->CellGeomCoeffValue(ii,jj,p1,p2);}
  //! Get the cell (ii,jj) geometric coefficient which is store in the 'position' place (i.e. powers and value)
  const double & CellGeomCoeffValue(const int & ii, const int & jj,
				    const int & position) {return Geom->CellGeomCoeffValue(ii,jj,position);}
  //@}

  //! @name Reconstruction order and number of rings
  //@{
  //! Get the number of rings of neighbour cells around the reconstructed cell.
  const int & Rings(void) const {return rings;}
  //! Get the number of rings of neighbour cells around the reconstructed cell used for smoothness indicator calculation.
  const int & RingsSI(void) const {return rings_SI;}
  //! Get the order of reconstruction set for this object.
  const int & RecOrder(void) const {return OrderOfReconstruction;}
  //@}

  //! @name Monotonicity info for high-order
  //@{
  //! Get the pointer to the double array of monotonicity flags.
  FlagType ** InadequateFit(void) { return LimitedCell;}
  //! Get the monotonicity flags for reconstruction of cell (ii,jj)
  const FlagType & CellInadequateFit(const int & ii, const int & jj) const { return LimitedCell[ii][jj];}
  //! Get the monotonicity flags for reconstruction of cell (ii,jj)
  FlagType & CellInadequateFit(const int & ii, const int & jj){ return LimitedCell[ii][jj];}
  /*!
   * Get the monotonicity flag for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
   * It is assumed that 'VarPosition' starts from ONE!
   */
  const int & CellInadequateFitValue(const int & ii, const int & jj,
				     const int VarPosition) const { return LimitedCell[ii][jj][VarPosition-1];}
  //! Get the monotonicity flag for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  int & CellInadequateFitValue(const int & ii, const int & jj,
			       const int VarPosition){ return LimitedCell[ii][jj][VarPosition-1];}
  const int & Previous_CellInadequateFitValue(const int & ii, const int & jj,
					      const int VarPosition) const { return PreviousLimitedCell[ii][jj][VarPosition-1];}
  int & Previous_CellInadequateFitValue(const int & ii, const int & jj,
					const int VarPosition){ return PreviousLimitedCell[ii][jj][VarPosition-1];}
  //@}

  //! @name Smoothness indicator
  //@{
  //! Get the pointer to the double array of reconstruction smoothness indicators.
  DoubleArrayType ** SmoothnessIndicator(void) { return SI;}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj)
  const DoubleArrayType & CellSmoothnessIndicator(const int & ii, const int & jj) const { return SI[ii][jj];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) 
  DoubleArrayType & CellSmoothnessIndicator(const int & ii, const int & jj){ return SI[ii][jj];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  const double & CellSmoothnessIndicatorValue(const int & ii, const int & jj,
					      const int & VarPosition) const { return SI[ii][jj][VarPosition-1];}
  //! Get the smoothness indicators for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  double & CellSmoothnessIndicatorValue(const int & ii, const int & jj,
					const int & VarPosition){ return SI[ii][jj][VarPosition-1];}
  //@}

  //! @name Pseudo-inverse of the LHS term in the CENO reconstruction
  //@{
  //! Get the pointer to the double array of reconstruction pseudo-inverse matrices.
  DenseMatrix ** LHS_Inv(void) {return CENO_LHS;}
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
  bool IsPseudoInversePreComputed(void) const { return _calculated_psinv; }
  //! Require update of the pseudo-inverse
  void MustUpdatePseudoInverse(void) { _calculated_psinv = false; }
  //! Return true if the pseudo-inverse related containers has been allocated.
  bool IsPseudoInverseAllocated(void) const { return _allocated_psinv; }
  //@}

  //! @name Geometric weights assigned to the cells that are part of the stencil
  //@{
  //! Get the pointer to the double array of geometric weights used for k-exact reconstruction.
  DoubleArrayType ** GeomWeights(void) {return CENO_Geometric_Weights;}
  //! Get the array of geometric weights for cell (ii,jj)
  DoubleArrayType & GeomWeights(const int & ii, const int & jj){return CENO_Geometric_Weights[ii][jj];}
  //! Get the array of geometric weights for cell (ii,jj)
  const DoubleArrayType & GeomWeights(const int & ii, const int & jj) const {return CENO_Geometric_Weights[ii][jj];}
  //! Get the geometric weight for the cell with the index 'CellPosition' which is part of the stencil of cell (ii,jj) 
  double & GeomWeightValue(const int & ii, const int & jj,
			   const int & CellPosition){return CENO_Geometric_Weights[ii][jj][CellPosition];}
  //! Get the geometric weight for the cell with the index 'CellPosition' which is part of the stencil of cell (ii,jj) 
  const double & GeomWeightValue(const int & ii, const int & jj,
				 const int & CellPosition) const {return CENO_Geometric_Weights[ii][jj][CellPosition];}
  //@}

  //! @name Block geometry
  //@{
  //! Get the associated block mesh
  const GeometryType* Geometry(void) const {return Geom;}
  //! Get the cell (ii,jj) of the associated block mesh
  GeometryType & CellGeometry(const int & ii, const int & jj){return Geom->Cell[ii][jj];}
  //! Get the cell (ii,jj) of the associated block mesh
  const GeometryType & CellGeometry(const int & ii, const int & jj) const {return Geom->Cell[ii][jj];}
  //! Get the centroid of cell (ii,jj) of the associated block mesh
  const Vector2D & CellCenter(const int & ii, const int & jj) const {return Geom->CellCentroid(ii,jj);}
  //! Get the X-coordinate of the (ii,jj) cell centroid of the associated block mesh
  const double & XCellCenter(const int & ii, const int & jj) const {return Geom->XCellCentroid(ii,jj);}
  //! Get the Y-coordinate of the (ii,jj) cell centroid of the associated block mesh
  const double & YCellCenter(const int & ii, const int & jj) const {return Geom->YCellCentroid(ii,jj);}
  //! Set the pointer to the associated geometry
  void SetGeometryPointer(GeometryType & Block){ Geom = &Block;}
  void AssociateGeometry(GeometryType & Block);
  //! Return true if constrained reconstruction is used anywhere in the block, otherwise return false
  const bool & IsConstrainedReconstructionRequired(void) const { return _constrained_block_reconstruction; }
  //! Return true if constrained reconstruction is used in the proximity of West block boundary, otherwise return false
  const bool IsWestConstrainedReconstructionRequired(void) const { return WestBnd.IsReconstructionConstrained(); }
  //! Return true if constrained reconstruction is used in the proximity of East block boundary, otherwise return false
  const bool IsEastConstrainedReconstructionRequired(void) const { return EastBnd.IsReconstructionConstrained(); }
  //! Return true if constrained reconstruction is used in the proximity of North block boundary, otherwise return false
  const bool IsNorthConstrainedReconstructionRequired(void) const { return NorthBnd.IsReconstructionConstrained(); }
  //! Return true if constrained reconstruction is used in the proximity of South block boundary, otherwise return false
  const bool IsSouthConstrainedReconstructionRequired(void) const { return SouthBnd.IsReconstructionConstrained(); }
  //! Get i-direction index of first interior cell
  const int & ICl_Grid(void) const {return Geom->ICl; }
  //! Get i-direction index of last interior cell
  const int & ICu_Grid(void) const {return Geom->ICu; }
  //! Get j-direction index of first interior cell
  const int & JCl_Grid(void) const {return Geom->JCl; }
  //! Get j-direction index of last interior cell
  const int & JCu_Grid(void) const {return Geom->JCu; }
  //! Get number of ghost cells for the associated block mesh
  const int & Nghost_Grid(void) const {return Geom->Nghost; }
  //@}

  //! @name Block boundary information
  //@{
  //! Get information for West boundary
  HighOrder2D_BlockBoundary & getWestBnd(void) { return WestBnd; }
  const HighOrder2D_BlockBoundary & getWestBnd(void) const { return WestBnd; }
  HighOrder2D_BlockBoundary & get_South_WestBnd(void) { return S_WestBnd; }
  const HighOrder2D_BlockBoundary & get_South_WestBnd(void) const { return S_WestBnd; }
  HighOrder2D_BlockBoundary & get_North_WestBnd(void) { return N_WestBnd; }
  const HighOrder2D_BlockBoundary & get_North_WestBnd(void) const { return N_WestBnd; }

  HighOrder2D_BlockBoundary & getEastBnd(void) { return EastBnd; }
  const HighOrder2D_BlockBoundary & getEastBnd(void) const { return EastBnd; }
  HighOrder2D_BlockBoundary & get_South_EastBnd(void) { return S_EastBnd; }
  const HighOrder2D_BlockBoundary & get_South_EastBnd(void) const { return S_EastBnd; }
  HighOrder2D_BlockBoundary & get_North_EastBnd(void) { return N_EastBnd; }
  const HighOrder2D_BlockBoundary & get_North_EastBnd(void) const { return N_EastBnd; }

  HighOrder2D_BlockBoundary & getNorthBnd(void) { return NorthBnd; }
  const HighOrder2D_BlockBoundary & getNorthBnd(void) const { return NorthBnd; }
  HighOrder2D_BlockBoundary & get_East_NorthBnd(void) { return E_NorthBnd; }
  const HighOrder2D_BlockBoundary & get_East_NorthBnd(void) const { return E_NorthBnd; }
  HighOrder2D_BlockBoundary & get_West_NorthBnd(void) { return W_NorthBnd; }
  const HighOrder2D_BlockBoundary & get_West_NorthBnd(void) const { return W_NorthBnd; }

  HighOrder2D_BlockBoundary & getSouthBnd(void) { return SouthBnd; }
  const HighOrder2D_BlockBoundary & getSouthBnd(void) const { return SouthBnd; }
  HighOrder2D_BlockBoundary & get_East_SouthBnd(void) { return E_SouthBnd; }
  const HighOrder2D_BlockBoundary & get_East_SouthBnd(void) const { return E_SouthBnd; }
  HighOrder2D_BlockBoundary & get_West_SouthBnd(void) { return W_SouthBnd; }
  const HighOrder2D_BlockBoundary & get_West_SouthBnd(void) const { return W_SouthBnd; }
  //@}

  //! @name Initialize container functions.
  //@{
  static int getMinimumNghost(const int &ReconstructionOrder); //!< return the minimum required number of ghost cells
  void SetRings(void);
  static int getStencilSize(const int &ReconstructionOrder);
  //! Return the stencil size for the current CENO reconstruction block.
  int getStencilSize(void) const { return getStencilSize(OrderOfReconstruction); }
  static int getTaylorDerivativesSize(const int &ReconstructionOrder);
  //! Return the number of Taylor derivatives for the current CENO reconstruction block.
  int getTaylorDerivativesSize(void) const { return getTaylorDerivativesSize(OrderOfReconstruction); }
  static int getNumberOfRings(int number_of_Taylor_derivatives);   //!< Return the required number of neighbour rings
  static int getNghostHighOrder(const int &ReconstructionOrder);
  //! Return the number of high-order ghost cells for the current CENO reconstruction block.
  int getNghostHighOrder(void) const { return getNghostHighOrder(OrderOfReconstruction); }
  const int & NghostHO(void) const { return Nghost_HO; }
  void ResetMonotonicityData(void);
  void ResetMonotonicityData(const int & ii, const int & jj);
  void ResetMonotonicityDataBackup(void);
  void ResetMonotonicityDataBackup(const int & ii, const int & jj);
  void InitializeMonotonicityVariables(const int & ii, const int & jj);
  void InitializeVariable(int ReconstructionOrder, GeometryType & Block,
			  const bool &_pseudo_inverse_allocation_ = false);
  void InitializeBasicVariable(int ReconstructionOrder, GeometryType & Block,
			       const bool &_pseudo_inverse_allocation_ = false);
  void SetReconstructionOrder(int ReconstructionOrder);
  void SetPropertiesHighOrderBlock(void);
  void CheckConsistencyOfGeometricSplineProperties(void);
  //@}

  //! @name Evaluate the polynomial interpolant.
  //@{
  //! Evaluate the interpolant at a given location (X_Coord,Y_Coord) for a specified solution variable (i.e. parameter),
  //  using the reconstruction of cell (ii,jj)
  double SolutionAtCoordinates(const int & ii, const int & jj, 
			       const double & X_Coord, const double & Y_Coord, const unsigned & parameter) const {
    return TD[ii][jj].ComputeSolutionFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj))[parameter];
  }
  //! Evaluate the interpolant at a given position vector for a specified solution variable (i.e. parameter),
  //  using the reconstruction of cell (ii,jj)
  double SolutionAtLocation(const int & ii, const int & jj,
			    const Vector2D &CalculationPoint, const unsigned & parameter) const {
    return SolutionAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y,parameter);
  }
  
  //! Evaluate the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj)
  Soln_State SolutionAtCoordinates(const int & ii, const int & jj,
				   const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeSolutionFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }

  //! Evaluate the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj)
  Soln_State SolutionStateAtCoordinates(const int & ii, const int & jj,
					const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeSolutionFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }

  //! Evaluate the interpolant at a given position vector for all solution variables, using the reconstruction of cell (ii,jj).
  Soln_State SolutionStateAtLocation(const int & ii, const int & jj, const Vector2D &CalculationPoint) const {
    return SolutionAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y);
  }
  
  //! Evaluate the X-gradient of the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State XGradientStateAtCoordinates(const int & ii, const int & jj,
					 const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeXGradientFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }
  //! Evaluate the X-gradient of the interpolant at a given position vector for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State XGradientStateAtLocation(const int & ii, const int & jj, const Vector2D &CalculationPoint) const {
    return XGradientStateAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y);
  }

  //! Evaluate the Y-gradient of the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State YGradientStateAtCoordinates(const int & ii, const int & jj,
					 const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeYGradientFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }
  //! Evaluate the Y-gradient of the interpolant at a given position vector for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State YGradientStateAtLocation(const int & ii, const int & jj, const Vector2D &CalculationPoint) const {
    return YGradientStateAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y);
  }

  //! Evaluate the entropy provided by the interpolant at a given location (X_Coord,Y_Coord),
  //  using the reconstruction of cell (ii,jj)
  double SolutionEntropyAtCoordinates(const int & ii, const int & jj, 
				      const double & X_Coord, const double & Y_Coord) const {
    return SolutionStateAtCoordinates(ii,jj,X_Coord,Y_Coord).s();
  }

  //! Evaluate the pressure provided by the interpolant at a given location (X_Coord,Y_Coord),
  //  using the reconstruction of cell (ii,jj)
  double SolutionPressureAtCoordinates(const int & ii, const int & jj, 
				      const double & X_Coord, const double & Y_Coord) const {
    return SolutionStateAtCoordinates(ii,jj,X_Coord,Y_Coord).pressure();
  }

  //! Evaluate the gradient in the normal direction of the interpolant 
  //  at given position coordinate for all solution variables, using the reconstruction of cell (ii,jj).
  Soln_State NormalGradientStateAtCoordinates(const int & ii, const int & jj,
					      const double & X_Coord, const double & Y_Coord,
					      const Vector2D &norm_dir) const;

  //! Evaluate the gradient in the normal direction of the interpolant 
  //  at given position coordinates for a specified solution variable (i.e. parameter), using the reconstruction of cell (ii,jj).
  double NormalGradientAtCoordinates(const int & ii, const int & jj,
				     const double & X_Coord, const double & Y_Coord,
				     const Vector2D &norm_dir, const unsigned & parameter) const;

  //! Evaluate the gradient in the normal direction of the interpolant 
  //  at a given position vector for all solution variables, using the reconstruction of cell (ii,jj).
  Soln_State NormalGradientStateAtLocation(const int & ii, const int & jj,
					   const Vector2D &CalculationPoint, const Vector2D &norm_dir) const {
    return NormalGradientStateAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y,norm_dir);
  }

  //! Evaluate the gradient in the normal direction of the interpolant 
  //  at a given position vector for a specified solution variable (i.e. parameter), using the reconstruction of cell (ii,jj).
  double NormalGradientAtLocation(const int & ii, const int & jj,
				  const Vector2D &CalculationPoint, const Vector2D &norm_dir,
				  const unsigned & parameter) const {
    return NormalGradientAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y,norm_dir,parameter);
  }
  //@}

  /*! @brief Integrate over the domain of the geometry associated with this high-order solution  */
  template<typename FO, class ReturnType>
  ReturnType IntegrateOverTheCell(const int &ii, const int &jj, const FO FuncObj,
				  const int & digits, ReturnType _dummy_param) const;

  
  /*! @brief Integrate the reconstructed polynomial of cell (ii,jj) over an arbitrary quad domain. */
  template<typename Node2DType>
  Soln_State IntegrateCellReconstructionOverQuadDomain(const int &ii, const int &jj,
						       const Node2DType &SW, const Node2DType &NW,
						       const Node2DType &NE, const Node2DType &SE) const;
  //! @name Reconstructions:
  //@{

  //! @name Block Level Reconstructions:
  //@{
  /*! @brief Compute the unlimited high-order solution reconstruction.  */
  template<class Soln_Block_Type>
  void ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk,
					      const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
												       const int &) const = 
					      &Soln_Block_Type::CellSolution);

  /*! @brief Compute the pseudo-inverse corresponding to the unlimited high-order solution reconstruction.  */
  void ComputeReconstructionPseudoInverse(void);

  /*! @brief Switch the previously detected non-smooth interpolants to limited piecewise linear reconstruction. */
  template<class Soln_Block_Type>
  void EnforceMonotonicityToNonSmoothInterpolants(Soln_Block_Type &SolnBlk,
						  const int &Limiter,
						  const Soln_State &
						  (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const = 
						  &Soln_Block_Type::CellSolution);

  /*! @brief Compute the high-order solution reconstruction (i.e. high-order, data analysis and monotonicity enforcement). */
  template<class Soln_Block_Type>
  void ComputeHighOrderSolutionReconstruction(Soln_Block_Type &SolnBlk,
					      const int &Limiter,
					      const Soln_State &
					      (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const = 
					      &Soln_Block_Type::CellSolution);
  //@} (Block Level Reconstructions)


  //! @name Cell Level Reconstructions:
  //@{
  //! @brief Return the type of reconstruction for a given cell (iCell,jCell)
  char getCellReconstructionType(const int& iCell, const int& jCell) const;
  char** const getReconstructionTypeMap(void) const {return ReconstructionTypeMap; }

  //! @brief Compute the unlimited high-order solution reconstruction of cell (iCell,jCell). 
  template<class Soln_Block_Type>
  void ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
					      const int &iCell, const int &jCell,
					      const bool & UseSpecialStencil,
					      const Soln_State & 
					      (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const  = 
					      &Soln_Block_Type::CellSolution );

  //! @brief Compute the unconstrained unlimited high-order solution reconstruction of cell (iCell,jCell). 
  template<class Soln_Block_Type>
  void ComputeUnconstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
							   const Soln_State & 
							   (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
							   const int &iCell, const int &jCell,
							   const IndexType & i_index, const IndexType & j_index);

  /*! @brief Compute the pseudo-inverse corresponding to the unlimited high-order solution reconstruction
    of cell (iCell,jCell).  */
  void ComputeCellReconstructionPseudoInverse(const int &iCell, const int &jCell,
					      const IndexType & i_index, const IndexType & j_index);
  /*! @brief Compute the pseudo-inverse or LHS matrix for cell (iCell,jCell) which is influenced
    by the presence of constrained boundaries.  */
  void ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(const int& iCell, const int& jCell);

  //! @brief Compute the piecewise linear solution reconstruction of cell (iCell,jCell). 
  template<class Soln_Block_Type>
  void ComputeLimitedPiecewiseLinearSolutionReconstruction(Soln_Block_Type &SolnBlk, 
							   const int &iCell, const int &jCell,
							   const int &Limiter,
							   const Soln_State &
							   (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const = 
							   &Soln_Block_Type::CellSolution);

  //! @brief Compute the constrained unlimited high-order solution reconstruction of cell (iCell,jCell). 
  template<class Soln_Block_Type>
  void ComputeConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
							 const Soln_State & 
							 (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
							 const int &iCell, const int &jCell,
							 const IndexType & i_index, const IndexType & j_index);

  //! @brief Compute the individually constrained unlimited high-order solution reconstruction of cell (iCell,jCell).
  template<class Soln_Block_Type>
  void ComputeIndividuallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
													   const int &) const,
								     const int &iCell, const int &jCell,
								     const IndexType & i_index, const IndexType & j_index,
								     const int & ConstrainedGQPs_West,
								     const int & ConstrainedGQPs_South,
								     const int & ConstrainedGQPs_East,
								     const int & ConstrainedGQPs_North);

  //! @brief Compute the relationally constrained unlimited high-order solution reconstruction of cell (iCell,jCell).
  template<class Soln_Block_Type>
  void ComputeRelationallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
								     const Soln_State & 
								     (Soln_Block_Type::*ReconstructedSoln)(const int &,
													   const int &) const,
								     const int &iCell, const int &jCell,
								     const IndexType & i_index, const IndexType & j_index,
								     const int & ConstrainedGQPs_West,
								     const int & ConstrainedGQPs_South,
								     const int & ConstrainedGQPs_East,
								     const int & ConstrainedGQPs_North);

  //! @brief Compute the relationally and individually constrained unlimited high-order solution reconstruction of cell (iCell,jCell).
  template<class Soln_Block_Type>
  void ComputeRelationallyAndIndividuallyConstrainedUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk, 
										    const Soln_State & 
										    (Soln_Block_Type::
										     *ReconstructedSoln)(const int &,
													 const int &) const,
										    const int &iCell, const int &jCell,
										    const IndexType & i_index,
										    const IndexType & j_index,
										    const int & ConstrainedGQPs_West,
										    const int & ConstrainedGQPs_South,
										    const int & ConstrainedGQPs_East,
										    const int & ConstrainedGQPs_North);

  //! @brief Compute the unlimited high-order reconstruction of those parameters that are unconstrained in a cell (iCell,jCell). 
  template<class Soln_Block_Type>
  void ComputeUnconstrainedUnlimitedSolutionReconstructionInConstrainedCell(Soln_Block_Type &SolnBlk, 
									    const Soln_State & 
									    (Soln_Block_Type::*ReconstructedSoln)(const int &,
														  const int &) const,
									    const int &iCell, const int &jCell,
									    const IndexType & i_index, const IndexType & j_index,
									    const IndexType & ParameterIndex);

  //! @brief Set the mean value conservation equations in the assemble matrix
  template<class Soln_Block_Type>
  void Set_MeanValueConservation_Equations(Soln_Block_Type & SolnBlk,
					   const Soln_State & 
					   (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const,
					   const int &iCell, const int &jCell,
					   const IndexType & i_index, const IndexType & j_index,
					   DenseMatrix & A, DenseMatrix & All_U,
					   const IndexType & ParameterIndex,
					   const int &RowConstraint,
					   const int &StartRow, const int &StartCol);

  //! @brief Set the mean value conservation equations in the LHS matrix of the k-exact reconstruction procedure
  void Set_LHS_MeanValueConservation_Equations(const int &iCell, const int &jCell,
					       const IndexType & i_index, const IndexType & j_index,
					       DenseMatrix & A, DoubleArrayType & GeometricWeights);

  //! @brief Set the individual constraint equations in the assemble matrix
  template<class Soln_Block_Type>
  void Generalized_IndividualConstraints_Equations(Soln_Block_Type & SolnBlk,
						   const int &iCell, const int &jCell,
						   Vector2DArray & Constraints_Loc,
						   Vector2DArray & Constraints_Normals,
						   BC_Type_Array & Constraints_BCs,
						   DenseMatrix & A, DenseMatrix & All_U,
						   const IndexType & ParameterIndex,
						   const int &StartRow, const int &StartCol);

  //! @brief Set the relational constraint equations in the assemble matrix
  template<class Soln_Block_Type>
  void Generalized_RelationalConstraints_Equations(Soln_Block_Type & SolnBlk,
						   const int &iCell, const int &jCell,
						   Vector2DArray & Constraints_Loc,
						   Vector2DArray & Constraints_Normals,
						   BC_Type_Array & Constraints_BCs,
						   const int & BC_Type,
						   DenseMatrix & A, DenseMatrix & All_U,
						   const IndexType & ParameterIndex,
						   const int &StartRow, const int &StartCol);

  //@} (Cell Level Reconstructions)

  //! @name Helper Functions:
  //@{
  /*! @brief Set the central stencil of cells used for reconstruction.  */
  void SetReconstructionStencil(const int &iCell, const int &jCell,
				IndexType & i_index, IndexType & j_index) const;
  void SetSpecialReconstructionStencil(const int &iCell, const int &jCell,
				       IndexType & i_index, IndexType & j_index) const;
  void SetConstrainedReconstructionStencil(const int &iCell, const int &jCell,
					   IndexType & i_index, IndexType & j_index) const;  
  void SetDeviatedReconstructionStencil(const int &iCell, const int &jCell,
					IndexType & i_index, IndexType & j_index,
					const int &rings,
					bool IsStencilExtended = true) const;
  void displayDeviatedReconstructionStencil(ostream & out,
					    const int &iCell, const int &jCell,
					    const int &rings) const;
  /*! @brief Set a central stencil of cells with a given extend for a particular cell. */
  void SetCentralStencil(const int &iCell, const int &jCell,
			 IndexType & i_index, IndexType & j_index,
			 const int &rings, int &StencilSize) const;
  /*! @brief Set a stencil with all cells used for the reconstructions of the cells in the passed index arrays. */
  void getEnlargedReconstructionStencil(const int &iCell, const int &jCell,
					IndexType & i_index, IndexType & j_index) const;
  //! @brief Gather the data on each cell edge for imposing the required constraints
  template<class Soln_Block_Type>
  void FetchDataConstraints(Soln_Block_Type & SolnBlk,
			    const int &iCell, const int &jCell,
			    int & BC_Type,
			    const int & parameter,
			    const int & ConstrainedGQPs_West,  const int * ConstraintBCs_W,
			    const int & ConstrainedGQPs_South, const int * ConstraintBCs_S,
			    const int & ConstrainedGQPs_East,  const int * ConstraintBCs_E,
			    const int & ConstrainedGQPs_North, const int * ConstraintBCs_N,
			    Vector2DArray & Constraints_Loc,
			    Vector2DArray & Constraints_Normals,
			    BC_Type_Array & Constraints_BCs);
  //@} (Helper Functions)

  //! @name CENO Analysis:
  //@{
  /*! @brief Compute smoothness indicator for all block cells. */
  template<class Soln_Block_Type>
  void ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk,
				  const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
											   const int &) const = 
				  &Soln_Block_Type::CellSolution);

  /*! @brief Compute smoothness data with a central stencil for a given cell (iCell,jCell). */
  template<class Soln_Block_Type>
  void ComputeCellSmoothnessDataWithCentralStencil(Soln_Block_Type &SolnBlk,
						   const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
													    const int &) const,
						   const int &iCell, const int &jCell);

  /*! @brief Compute smoothness data with a deviated stencil for a given cell (iCell,jCell). */
  template<class Soln_Block_Type>
  void ComputeCellSmoothnessDataWithDeviatedStencil(Soln_Block_Type &SolnBlk,
						    const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
													     const int &) const,
						    const int &iCell, const int &jCell);

  /*! @brief Compute smoothness data for a given cell (iCell,jCell) which is part of a constrained block. */
  template<class Soln_Block_Type>
  void ComputeCellSmoothnessDataForConstrainedBlock(Soln_Block_Type &SolnBlk,
						    const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
													     const int &) const,
						    const int &iCell, const int &jCell);
  
  /*! @brief Compute smoothness indicator for a specific cell (iCell,jCell). */
  template<class Soln_Block_Type>
  void ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk,
				  const Soln_State & (Soln_Block_Type::*ReconstructedSoln)(const int &,
											   const int &) const,
				  const int &iCell, const int &jCell,
				  const IndexType & i_index, const IndexType & j_index, const int & StencilSize);
  //! @brief Check whether any of the high-order interpolants of (iCell,jCell) cell is flagged as non-smooth.
  const bool IsThereAnyNonSmoothHighOrderReconstruction(const int &iCell, const int &jCell) const;
  //@}

  //! @name Access to the first-order derivatives from memory pool:
  //@{
  //! Get dUdx
  const Soln_State & get_dUdx(void) { return dUdx; }
  //! Get dUdy
  const Soln_State & get_dUdy(void) { return dUdy; }
  //@}

  //! @name Functions for positivity analysis of elliptic discretization:
  //@{
  //! @brief Compute Green-Gauss gradient at inter-cellular face for first parameter
  template<class Soln_Block_Type>
  void GreenGauss_FaceGradient_CentroidPathCartesianMesh(Soln_Block_Type &SolnBlk,
							 const int &iCell, const int &jCell,
							 const int &Face,
							 Vector2D & GradU_face,
							 const Soln_State &
							 (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const = 
							 &Soln_Block_Type::CellSolution);
  //@}

  //! @name Error Evaluation:
  //@{

  //! @name Access to the error data:
  //@{
  const double & L1(void) const {return ErrorL1; }       //!< return the L1 error
  const double & L2(void) const {return ErrorL2; }       //!< return the L2 error
  const double & LMax(void) const {return ErrorMax; }    //!< return the LMax error
  const double & BlockArea(void) const {return TotalBlockArea; } //!< return the area of the block used to calculate errors
  
  double BlockL1Norm(void) { return ErrorL1/TotalBlockArea; }	      //!< return the L1 error norm for the block
  double BlockL2Norm(void) { return sqrt(ErrorL2/TotalBlockArea); }   //!< return the L2 error norm for the block
  double BlockLMaxNorm(void) { return ErrorMax; }                     //!< return the LMax error norm for the block
  
  const unsigned int & UsedCells(void) const {return CellsUsed;}  //!< return the number of used cells for error calculation
  //@}

  //! @name Block Level Error Computation:
  //@{
  /*! @brief Compute the integral of the solution error function (i.e. L1, L2 and LMax errors) for the whole block */
  template<typename Function_Object_Type>
  void ComputeSolutionErrors(const Function_Object_Type FuncObj, const unsigned &parameter,
			     const int & digits = 8);

  /*! @brief Compute the integral of the error of two high-order reconstructions (i.e. L1, L2 and LMax norm) for the whole block */
  void ComputeReconstructionsErrors(const HighOrder2D<Soln_State> & Obj, const unsigned &parameter,
				    const int & digits = 8);
  //@} (Block Level Errors)

  //! @name Cell Level Errors Computation:
  //@{
  /*! @brief Compute the integral of the solution error function (i.e. L1 norm) for cell (iCell,jCell) */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL1(const int &iCell, const int &jCell,
				const Function_Object_Type FuncObj, const unsigned &parameter,
				const int & digits = 8);

  /*! @brief Compute the integral of the solution entropy error function (i.e. L1 norm) for cell (iCell,jCell) */
  template<typename Function_Object_Type>
  double ComputeSolutionEntropyErrorL1(const int &iCell, const int &jCell,
				       const Function_Object_Type FuncObj,
				       const int & digits = 8);

  /*! @brief Compute the integral of the error of two high-order reconstructions (i.e. L1 norm) for cell (iCell,jCell) */
  double ComputeReconstructionsErrorL1(const int &iCell, const int &jCell,
				       const HighOrder2D<Soln_State> & Obj, const unsigned &parameter,
				       const int & digits = 8);

  /*! @brief Compute the integral of the squared solution error function (i.e. L2 norm) for cell (iCell,jCell) */
  template<typename Function_Object_Type>
  double ComputeSolutionErrorL2(const int &iCell, const int &jCell,
				Function_Object_Type FuncObj, const unsigned &parameter,
				const int & digits = 8);

  /*! @brief Compute the integral of the squared solution entropy error function (i.e. L2 norm) for cell (iCell,jCell) */
  template<typename Function_Object_Type>
  double ComputeSolutionEntropyErrorL2(const int &iCell, const int &jCell,
				       Function_Object_Type FuncObj,
				       const int & digits = 8);

  /*! @brief Compute the integral of the squared error of two high-order reconstructions (i.e. L1 norm) for cell (iCell,jCell) */
  double ComputeReconstructionsErrorL2(const int &iCell, const int &jCell,
				       const HighOrder2D<Soln_State> & Obj, const unsigned &parameter,
				       const int & digits = 8);
  //@} (Cell Level Errors)

  //@} (Error Evaluation)

  //! @name AMR functions:
  //@{
  template<class Soln_Block_Type>
  double AMR_Criteria_Based_On_Minimum_Smoothness_Indicator(Soln_Block_Type &SolnBlk);
  //@}

  //! @name Input/Output functions:
  //@{
  void Output_Object(ostream & out_file) const;
  void Read_Object(istream & in_file);
  void outputReconstructionTypeMap(ostream & out_file) const;
  //@}

  //! @name Broadcast functions (MPI):
  //@{
  void Broadcast_HighOrder_Data(GeometryType & Block_Geometry);
#ifdef _MPI_VERSION
  void Broadcast_HighOrder_Data(MPI::Intracomm &Communicator, 
				const int &Source_CPU,
				GeometryType & Block_Geometry);
#endif
  //@}

  //! @name Friend functions:
  //@{
  friend ostream & operator<< <Soln_State> (ostream & os, const HighOrder2D<Soln_State> & Obj);
  friend istream & operator>> <Soln_State> (istream & os, HighOrder2D<Soln_State> & Obj);
  //@}

protected:
  
private:
  //! @name Internal indexes:
  //@{
  int Ni;		       //!< Number of high-order objects in i-direction
  int Nj;		       //!< Number of high-order objects in j-direction
  int Ng; 		       //!< Number of block ghost cells 
  int ICl,		       //!< Index of first interior cell in i-direction 
    ICu,		       //!< Index of last interior cell in i-direction 
    JCl,		       //!< Index of first interior cell in j-direction 
    JCu;                       //!< Index of last interior cell in j-direction 
  int OrderOfReconstruction;   //!< The order of reconstruction of the high-order object.
  int Nghost_HO;	       //!< Number of ghost cells in which high-order reconstruction is performed. 
  int StartI,		       //!< Index of the first cell in i-direction in which NO constrained reconstruction is performed. 
    EndI,                      //!< Index of the last cell in i-direction in which NO constrained reconstruction is performed. 
    StartJ,                    //!< Index of the first cell in j-direction in which NO constrained reconstruction is performed. 
    EndJ;                      //!< Index of the last cell in j-direction in which NO constrained reconstruction is performed. 
  //@}

  //! @name Memory allocation flags:
  //@{
  bool _allocated_block;       //!< Flag indicating if the containers at block level has been allocated or not.
  bool _allocated_cells;       //!< Flag indicating if the containers at cell level has been allocated or not.
  bool _allocated_psinv;       //!< Flag indicating if the pseudo-inverse related containers have been allocated or not. 
  //@}

  //! @name Member functions for limiter freezing:
  //@{
  bool _freeze_limiter;	       //!< Flag indicating if the limiter value must be frozen. Set based on the solution block flag.
  //@}

  //! @name Reconstruction containers:
  //@{
  DerivativesContainer **TD;   //!< High-order TaylorDerivatives
  DoubleArrayType **SI;        //!< The values of the smoothness indicator calculated for each reconstructed variable
  FlagType **LimitedCell;      //!< Monotonicity flag: Values --> OFF - high-order reconstruction,
                               //                                  ON - limited linear reconstruction
  FlagType **PreviousLimitedCell;      //!< Copy of the LimitedCell variable from the previous reconstruction

  DenseMatrix **CENO_LHS;      //!< Storage for the pseudo-inverse of the LHS term in the CENO reconstruction.
  DoubleArrayType **CENO_Geometric_Weights;   //!< Storage for the geometric weights used in CENO reconstruction.
  /*!
   * Storage for the reconstruction type of cells that are part of a block that is flagged as constrained.
   * For unconstrained blocks this variable is NULL!
   * Possible types: 
   *    n - no reconstruction (i.e. cells for which no reconstruction is carried out)
   *    r - regular reconstruction (i.e. cells that use the central stencil for reconstruction)
   *    m - modified reconstruction (i.e. cells that have the modified reconstruction stencil and no constraints)
   *    c - constrained reconstruction (i.e. cells that have the modified reconstruction stencil and constraints)
   */
  char **ReconstructionTypeMap;
  //@}

  //! @name Other internal variables/flags:
  //@{
  int rings;                   //!< Number of rings used to generate the reconstruction stencil
  int rings_SI;		       //!< Number of rings used to compute the smoothness indicator
  bool _calculated_psinv;      //!< Flag to indicate whether the pseudo-inverse has been already computed or not.
  int _si_calculation;         /*!< Flag to store the method for smoothness indicator calculation. Set to the 
				* same value as CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS. */
  //@}

  //! @name Correspondent grid trackers:
  //@{
  unsigned int ObserverInteriorCellGeometryState; //!< Observer to memorise the state of the interior grid.
  unsigned int ObserverGhostCellGeometryState;	  //!< Observer to memorise the state of the ghost cell layers.
  unsigned int ObserverCornerGhostCellGeometryState; //!< Observer to memorize the state of the corner ghost cells.
  //@}

  //! @name Internal information in conjunction with the geometry:
  //@{
  GeometryType* Geom;          //!< Pointer to the associated geometry which is a 2D mesh block
  bool _constrained_block_reconstruction; /*!< Flag to indicate whether constrained reconstruction is carried out for
					   * this block on any of the block boundaries. */
  HighOrder2D_BlockBoundary WestBnd,	  //!< Storage for information related to the West high-order block boundary
    S_WestBnd,                  //!< Storage for information related to the South extension of West high-order block boundary
    N_WestBnd,                  //!< Storage for information related to the North extension of West high-order block boundary
    EastBnd,			//!< Storage for information related to the East high-order block boundary
    S_EastBnd,			//!< Storage for information related to the South extension of East high-order block boundary
    N_EastBnd,			//!< Storage for information related to the North extension of East high-order block boundary
    NorthBnd,                   //!< Storage for information related to the North high-order block boundary
    W_NorthBnd,                 //!< Storage for information related to the West extension of North high-order block boundary
    E_NorthBnd,                 //!< Storage for information related to the East extension of North high-order block boundary
    SouthBnd,                   //!< Storage for information related to the South high-order block boundary
    W_SouthBnd,                 //!< Storage for information related to the West extension of South high-order block boundary
    E_SouthBnd;                 //!< Storage for information related to the East extension of South high-order block boundary
  //! Reset block boundaries information
  void ResetBlockBoundaries(void);
  //@}

  //! Get the number of variables in the solution state
  int NumberOfVariables(void) const {return Soln_State::NumVar(); }

  //! Allocate memory at the cell level based on the order of reconstruction
  void allocate_CellMemory(const int &ReconstructionOrder, const bool &_pseudo_inverse_allocation_);
  //! Allocate memory for the reconstruction type map
  void allocate_ReconstructionTypeMap(void);
  //! @brief Build reconstruction type map
  void BuildReconstructionTypeMap(void);

  //! @name Reconstruction memory pools:
  //@{
  // === Helper variables (i.e. memory pools which are overwritten for each cell in the reconstruction process) ===
  Vector2DArray DeltaCellCenters;    // Storage for distances between centroids.
  IndexType i_index, j_index;	     // Storage for indexes of cells that are part of the reconstruction stencil.
  DenseMatrix A;		     // Storage for the LHS matrix of the k-Exact reconstruction.
  DoubleArrayType GeometricWeights;  // Storage for the geometric weights calculated in the k-Exact reconstruction.
  DenseMatrix All_Delta_U;	     // Storage for the RHS matrix (i.e. solution dependent) of the k-Exact reconstruction.
  ColumnVector Delta_U;              // Storage for a particular column of the RHS matrix.
  ColumnVector X;                    // Storage for the solution to the least-square problem.

  // === Helper variables for the limited piecewise linear solution reconstruction ===
  Soln_State U_ave, dUdx, dUdy;
  int I_Index[8], J_Index[8];
  double geom_weights[8];
  Vector2D dX[8];

  // === Helper variables for constrained reconstruction
  // (i.e. memory pools which are overwritten for each cell in the constrained reconstruction process) ===
  IndexType i_index_ave, j_index_ave;	/* Storage for indexes of cells that contribute to the reconstruction stencil 
					   with conservation of average solution values */
  Vector2DArray Constraints_Loc;	// Storage for the locations where constraints are imposed
  Vector2DArray Constraints_Normals;	// Storage for the normal vectors at the locations where constraints are imposed
  BC_Type_Array Constraints_BCs;	// Storage for the high-order boundary conditions that are imposed as constraints

  Vector2DArray Approx_Constraints_Loc;	     // Storage for the locations where approximate constraints are imposed
  Vector2DArray Approx_Constraints_Normals;  // Storage for the normal vectors at the locations where approx. constraints are imposed
  BC_Type_Array Approx_Constraints_BCs;	     // Storage for the high-order boundary conditions that are imposed as approx. constraints

  DenseMatrix A_Assembled;	 // Storage for the LHS matrix of the constrained k-Exact reconstruction.
  DenseMatrix All_U_Assembled;   // Storage for the RHS matrix (i.e. solution dependent) of the constrained k-Exact reconstruction.
  DenseMatrix X_Assembled;       // Storage for the solution to the constrained least-square problem.
  //@}

  //! @name Error calculation variables and internal functions:
  //@{
  double ErrorL1, ErrorL2, ErrorMax; //!< errors/error norms
  double TotalBlockArea;             //!< the area of the block used for error calculation
  unsigned int CellsUsed;	     //!< the number of cells used for accuracy assessment
  //! Set all error variables to zero
  void ResetErrors(void){ ErrorL1 = 0.0; ErrorL2 = 0.0; ErrorMax = 0.0; TotalBlockArea = 0.0; CellsUsed = 0;}
  //@}

  //! @name Functions related to the memory pool piecewise linear solution interpolant:
  //@{
  //! Evaluate the linear interpolant at a given location (X_Coord,Y_Coord) 
  //  for a specified solution variable (i.e. parameter) and a cell (ii,jj)
  double UnlimitedLinearSolutionAtCoordinates(const int & ii, const int & jj, 
					      const double & X_Coord, const double & Y_Coord, const unsigned & parameter) const {
    return U_ave[parameter] + dUdx[parameter]*(X_Coord - XCellCenter(ii,jj)) + dUdy[parameter]*(Y_Coord - YCellCenter(ii,jj));
  }
  //! Evaluate the linear interpolant at a given position vector for a 
  // specified solution variable (i.e. parameter) and a cell (ii,jj)
  double UnlimitedLinearSolutionAtLocation(const int & ii, const int & jj,
					   const Vector2D &CalculationPoint, const unsigned & parameter) const {
    return UnlimitedLinearSolutionAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y,parameter);
  }
  //! Evaluate the limiter
  double CalculateLimiter(double *uQuad, const int &nQuad,
			  const double &u0, const double &u0Min, const double &u0Max, const int &Limiter);
  //@}

  //! @name Functions and variables related to calculation of the smoothness indicator
  //@{
  int StencilSize_SmoothnessIndicator;
  //! Mean solution of the cell for which the smoothness indicator is computed. It's used as a reference.
  Soln_State MeanSolution;
  //! Sum of the squares of differences between the reconstructions of stencil cells evaluated at their centroids and MeanSolution.
  Soln_State SS_Regression;
  /*! Sum of the squares of difference between the reconstruction of cell (iCell,jCell) and 
    the reconstructions of the stencil cells evaluated at their centroids. */
  Soln_State SS_Residual;
  //! Intermediate variables.
  Soln_State Max_SS_Regression, Temp_Regression, Temp_Residual, Min_Mean_Solution;
  //! Normalization state characteristic to each solution state and cell
  Soln_State NormalizationState;
  //! Coefficient used to compute the smoothness indicator. It is based on the degrees of freedom (i.e. number of derivatives) and the size of the reconstruction stencil
  double AdjustmentCoeff;
  //! Local variable
  Vector2D DeltaCentroids;
  //! Local variable
  double GeomWeightSI;

  //! @brief Flag (iCell,jCell) cell with non-smooth high-order interpolants.
  void FlagCellReconstructionsAsNonSmooth(const int &iCell, const int &jCell);
  /*! @brief Check the smoothness condition and flag the non-smooth interpolants for a given cell.*/
  void AssessInterpolantsSmoothness(const int &iCell, const int &jCell);
  //@}

};

/******************************************************
 * Initialize the HighOrder2D class static variables  *
 *****************************************************/


/****************************************************
 * Implement the HighOrder2D class member functions *
 ***************************************************/
//! Default Constructor
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(void):
  Ni(0), Nj(0), Ng(0),
  ICl(0), ICu(0), JCl(0), JCu(0),
  OrderOfReconstruction(-1), Nghost_HO(0),
  StartI(0), EndI(0), StartJ(0), EndJ(0),
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  _freeze_limiter(false),
  TD(NULL), SI(NULL), LimitedCell(NULL), PreviousLimitedCell(NULL),
  rings(0), rings_SI(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL),
  ReconstructionTypeMap(NULL),
  Geom(NULL), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS),
  _constrained_block_reconstruction(false),
  AdjustmentCoeff(0.0),
  ObserverInteriorCellGeometryState(0), ObserverGhostCellGeometryState(0), ObserverCornerGhostCellGeometryState(0)
{
  //
}

//! Main Constructor
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(int ReconstructionOrder,
				     GeometryType & Block,
				     const bool &_pseudo_inverse_allocation_):
  Ni(0), Nj(0), Ng(0),
  ICl(0), ICu(0), JCl(0), JCu(0),
  OrderOfReconstruction(-1), Nghost_HO(0),
  StartI(0), EndI(0), StartJ(0), EndJ(0),
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  _freeze_limiter(false),
  TD(NULL), SI(NULL), LimitedCell(NULL), PreviousLimitedCell(NULL),
  rings(0), rings_SI(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL), 
  ReconstructionTypeMap(NULL),
  Geom(&Block), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS),
  _constrained_block_reconstruction(false),
  AdjustmentCoeff(0.0),
  ObserverInteriorCellGeometryState(Block.getInteriorStateTracker()),
  ObserverGhostCellGeometryState(Block.getGhostStateTracker()),
  ObserverCornerGhostCellGeometryState(Block.getCornerGhostStateTracker())
{
  // Use the grid to get the number of interior block cells and the number of available ghost cells
  allocate(ICu_Grid()-ICl_Grid()+1,
	   JCu_Grid()-JCl_Grid()+1,
	   Nghost_Grid(),
	   _pseudo_inverse_allocation_,
	   ReconstructionOrder);

  // Set properties of block boundaries (i.e. which boundaries are constrained/unconstrained and which are opaque/transparent)
  SetPropertiesHighOrderBlock();

  // Build the reconstruction type map
  BuildReconstructionTypeMap();

}

//! Copy constructor 
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(const HighOrder2D<SOLN_STATE> & rhs)
  : Ni(0), Nj(0), Ng(0),
    ICl(0), ICu(0), JCl(0), JCu(0),
    OrderOfReconstruction(-1), Nghost_HO(0),
    _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
    _freeze_limiter(false),
    TD(NULL), SI(NULL), LimitedCell(NULL), PreviousLimitedCell(NULL),
    rings(0), rings_SI(0), _calculated_psinv(false),
    CENO_LHS(NULL), CENO_Geometric_Weights(NULL),
    ReconstructionTypeMap(NULL),
    Geom(rhs.Geom), _si_calculation(rhs._si_calculation),
    AdjustmentCoeff(rhs.AdjustmentCoeff),
    ObserverInteriorCellGeometryState(rhs.ObserverInteriorCellGeometryState),
    ObserverGhostCellGeometryState(rhs.ObserverGhostCellGeometryState),
    ObserverCornerGhostCellGeometryState(rhs.ObserverCornerGhostCellGeometryState)
{

  int i,j;

  // check if the rhs has block memory allocated
  if (rhs._allocated_block){

    // allocate memory for the new container
    allocate(rhs.Ni - 2*rhs.Ng,
	     rhs.Nj - 2*rhs.Ng,
	     rhs.Ng,
	     rhs._allocated_psinv,
	     rhs.OrderOfReconstruction);
    
    // check if the rhs has cell memory allocated
    if (rhs._allocated_cells){

      // copy the cell containers
      for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
	for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) {
	  TD[i][j] = rhs.TD[i][j];
	  SI[i][j] = rhs.SI[i][j];
	  LimitedCell[i][j] = rhs.LimitedCell[i][j];
	  PreviousLimitedCell[i][j] = rhs.PreviousLimitedCell[i][j];

	  // copy the pseudo-inverse data
	  if (rhs._allocated_psinv){
	    CENO_LHS[i][j] = rhs.CENO_LHS[i][j];
	    CENO_Geometric_Weights[i][j] = rhs.CENO_Geometric_Weights[i][j];
	  } // endif (rhs._allocated_psinv)
	}//endfor
      }//endfor

    } else {
      deallocate_CellMemory();
    } //endif (rhs._allocated_cells)

    // set range of cells without constrained reconstruction
    StartI = rhs.StartI;
    EndI = rhs.EndI;
    StartJ = rhs.StartJ;
    EndJ = rhs.EndJ;

    // set the type of reconstruction required by the geometry setup
    _constrained_block_reconstruction = rhs._constrained_block_reconstruction;

    // copy properties of block boundaries
    WestBnd = rhs.WestBnd;
    S_WestBnd = rhs.S_WestBnd;
    N_WestBnd = rhs.N_WestBnd;

    EastBnd = rhs.EastBnd;
    S_EastBnd = rhs.S_EastBnd;
    N_EastBnd = rhs.N_EastBnd;
    
    NorthBnd = rhs.NorthBnd;
    E_NorthBnd = rhs.E_NorthBnd;
    W_NorthBnd = rhs.W_NorthBnd;

    SouthBnd = rhs.SouthBnd;
    E_SouthBnd = rhs.E_SouthBnd;
    W_SouthBnd = rhs.W_SouthBnd;

    // determine the reconstruction type of each cell for a constrained block
    BuildReconstructionTypeMap();
    
  }//endif (rhs._allocated_block)
}

//! Assignment operator
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE> & HighOrder2D<SOLN_STATE>::operator=(const HighOrder2D<SOLN_STATE> & rhs){

  int i,j;

  // Handle self-assignment:
  if (this == & rhs) return *this;

  // check if the rhs has block memory allocated
  if (rhs._allocated_block){

    // allocate memory for the new container
    allocate(rhs.Ni - 2*rhs.Ng,
	     rhs.Nj - 2*rhs.Ng,
	     rhs.Ng,
	     rhs._allocated_psinv,
	     rhs.OrderOfReconstruction);

    // make sure that the two objects have the same execution mode
    if (_si_calculation != rhs._si_calculation){
      deallocate();
      throw runtime_error("HighOrder2D<SOLN_STATE>::operator=() ERROR! The object cannot be assigned due to incompatibilities between the CENO_Execution_Mode class settings and the object settings");
    }

    // set the smoothness indicator adjustment coefficient
    AdjustmentCoeff = rhs.AdjustmentCoeff;
    
    // check if the rhs has cell memory allocated
    if (rhs._allocated_cells){

      // copy the cell containers
      for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
	for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) {
	  TD[i][j] = rhs.TD[i][j];
	  SI[i][j] = rhs.SI[i][j];
	  LimitedCell[i][j] = rhs.LimitedCell[i][j];
	  PreviousLimitedCell[i][j] = rhs.PreviousLimitedCell[i][j];

	  // copy the pseudo-inverse data
	  if (rhs._allocated_psinv){
	    CENO_LHS[i][j] = rhs.CENO_LHS[i][j];
	    CENO_Geometric_Weights[i][j] = rhs.CENO_Geometric_Weights[i][j];
	  } // endif (rhs._allocated_psinv)
	}//endfor
      }//endfor

    } else {
      deallocate_CellMemory();
    } //endif (rhs._allocated_cells)

    // set range of cells without constrained reconstruction
    StartI = rhs.StartI;
    EndI = rhs.EndI;
    StartJ = rhs.StartJ;
    EndJ = rhs.EndJ;

    // set the type of reconstruction required by the geometry setup
    _constrained_block_reconstruction = rhs._constrained_block_reconstruction;

    // copy properties of block boundaries
    WestBnd = rhs.WestBnd;
    S_WestBnd = rhs.S_WestBnd;
    N_WestBnd = rhs.N_WestBnd;

    EastBnd = rhs.EastBnd;
    S_EastBnd = rhs.S_EastBnd;
    N_EastBnd = rhs.N_EastBnd;
    
    NorthBnd = rhs.NorthBnd;
    E_NorthBnd = rhs.E_NorthBnd;
    W_NorthBnd = rhs.W_NorthBnd;

    SouthBnd = rhs.SouthBnd;
    E_SouthBnd = rhs.E_SouthBnd;
    W_SouthBnd = rhs.W_SouthBnd;

    // determine the reconstruction type of each cell for a constrained block
    BuildReconstructionTypeMap();

  }//endif (rhs._allocated_block)

  // Set grid observers
  ObserverInteriorCellGeometryState = rhs.ObserverInteriorCellGeometryState;
  ObserverGhostCellGeometryState = rhs.ObserverGhostCellGeometryState;
  ObserverCornerGhostCellGeometryState = rhs.ObserverCornerGhostCellGeometryState;

  return *this;
}

// allocate()
/*! Allocate memory for the high-order object.
 *
 * \param NC_Idir number of cells in i-direction
 * \param NC_Jdir number of cells in j-direction
 * \param Nghost  number of block ghost cells
 * \param ReconstructionOrder the order of reconstruction for this high-order object.
 *        If this number is not specified the implicit value is -1, which corresponds to no memory allocation.
 * 
 * \todo Improve the memory allocation based on the runtime situations and reconstruction order.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::allocate(const int &NC_IDir,
				       const int &NC_JDir,
				       const int &Nghost,
				       const bool &_pseudo_inverse_allocation_,
				       int ReconstructionOrder){

  int i,j;

  // Check conditions
  if (NC_IDir < 2 || NC_JDir < 2 || Nghost < 1 || ReconstructionOrder < -1){
    throw runtime_error("HighOrder2D<SOLN_STATE>::allocate() ERROR! Inconsistent dimensions!");
  }

  // Don't allocate any memory and deallocate any previously allocated one if reconstruction order is -1
  if (ReconstructionOrder == -1){
    deallocate();
    return;
  }

  if (CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING){
    // This method is not implemented yet but it's left as a possible development stream.
    // In the current implementation the 'CENO_Execution_Mode' class should refuse to set this flag to ON.

  } else {

    // Check if the new block dimensions are different than the currently allocated ones
    if ( Ni != (NC_IDir+2*Nghost) || Nj != (NC_JDir+2*Nghost) || Ng != Nghost ||
	 _si_calculation != CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS ){

      // free the memory if there is memory allocated
      deallocate();

      /* check consistency relationship between Nghost and 
	 the minimum number of ghost cells required for the
	 provided ReconstructionOrder */
      if ( Nghost < getMinimumNghost(ReconstructionOrder) ){
	throw runtime_error("HighOrder2D<SOLN_STATE>::allocate() ERROR! Too few ghost cells provided for the required reconstruction!");
      };

      // allocate new memory 
      Ni = NC_IDir+2*Nghost; Nj = NC_JDir+2*Nghost;
      Ng = Nghost;
      ICl = Ng; ICu = Ni-Ng-1;
      JCl = Ng; JCu = Nj-Ng-1;
    
      // Allocate memory at block level
      TD = new DerivativesContainer* [Ni];
      SI = new DoubleArrayType* [Ni];
      LimitedCell = new FlagType* [Ni];
      PreviousLimitedCell = new FlagType* [Ni];
      if (_pseudo_inverse_allocation_){
	CENO_LHS = new DenseMatrix* [Ni];
	CENO_Geometric_Weights = new DoubleArrayType* [Ni];
      }

      for (i = 0; i < Ni ; ++i){
	TD[i] = new DerivativesContainer [Nj];
	SI[i] = new DoubleArrayType [Nj];
	LimitedCell[i] = new FlagType [Nj];
	PreviousLimitedCell[i] = new FlagType [Nj];
	if (_pseudo_inverse_allocation_){
	  CENO_LHS[i] = new DenseMatrix [Nj];
	  CENO_Geometric_Weights[i] = new DoubleArrayType [Nj];
	}
      }// endfor

      // Block memory allocated
      _allocated_block = true;

      // Allocate memory and initialize at cell level.
      allocate_CellMemory(ReconstructionOrder,
			  _pseudo_inverse_allocation_);

    } else if ( ReconstructionOrder != OrderOfReconstruction ) { 

      // re-allocate only the cell memory
      allocate_CellMemory(ReconstructionOrder,
			  _pseudo_inverse_allocation_);

    } else {
      // Reset the backup of the monotonicity data
      ResetMonotonicityDataBackup();

    }//endif

  }//endif
}

// allocate_CellMemory()
/*!
 * Allocate memory for all cells and set the variables
 * that depend on the order of reconstruction.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::allocate_CellMemory(const int &ReconstructionOrder, const bool &_pseudo_inverse_allocation_){

  int i,j;
  int StencilSize;

  // Set the new reconstruction order
  OrderOfReconstruction = ReconstructionOrder;

  // Set the Nghost_HO based on the OrderOfReconstruction
  Nghost_HO = getNghostHighOrder();

  // Set the number of neighbour rings based on the OrderOfReconstruction
  SetRings();

  // Store stencil size
  StencilSize = getStencilSize();

  // Allocate memory and initialize containers at cell level.
  for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
    for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) {    
	
      // Set Taylor derivatives
      TD[i][j].GenerateContainer(OrderOfReconstruction);

      // Set smoothness indicator and monotonicity flag
      InitializeMonotonicityVariables(i,j);

      // Allocate pseudo-inverse data
      // Note: There is no need to initialize these containers here!
      if (_pseudo_inverse_allocation_){
	CENO_LHS[i][j].newsize(StencilSize - 1, TD[i][j].size() - 1);
	CENO_Geometric_Weights[i][j].assign(StencilSize, 0.0);
      }

    } /* endfor */
  }/* endfor */

  // Set the reconstruction helper variables
  DeltaCellCenters.assign(StencilSize, 0.0); // This variable is overwritten for each cell
  i_index.assign(StencilSize, 0); // This variable is overwritten for each cell
  j_index.assign(StencilSize, 0); // This variable is overwritten for each cell
  A.newsize(StencilSize - 1, NumberOfTaylorDerivatives() - 1);
  GeometricWeights.assign(StencilSize, 0.0);
  All_Delta_U.newsize(StencilSize - 1, NumberOfVariables());
  Delta_U.newsize(StencilSize - 1);
  X.newsize(NumberOfTaylorDerivatives() - 1);

  // Set the constrained reconstruction helper variables
  i_index_ave.reserve(StencilSize);
  j_index_ave.reserve(StencilSize);
  Constraints_Loc.reserve(12);
  Constraints_Normals.reserve(12);
  Constraints_BCs.reserve(12);
  Approx_Constraints_Loc.reserve(50);
  Approx_Constraints_Normals.reserve(50);
  Approx_Constraints_BCs.reserve(50);
  X_Assembled.newsize(NumberOfTaylorDerivatives(), NumberOfVariables());

  // Confirm allocation
  _allocated_cells = true;
  if (_pseudo_inverse_allocation_){
    _allocated_psinv = true;
    _calculated_psinv = false;
  }

  // Remember the smoothness indicator calculation method which was used for the current setup.
  _si_calculation = CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS;
}

/*!
 * Allocate memory for the reconstruction type 
 * of each cell only if the block is flagged as constrained.
 * If the block is not constrained and memory has 
 * been previously allocated to ReconstructionTypeMap variable
 * this will be deallocated.
 *
 * \note Call this routine AFTER the block type has been decided!
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::allocate_ReconstructionTypeMap(void){

  if (_allocated_block && IsConstrainedReconstructionRequired()){
    if (ReconstructionTypeMap == NULL){
      // Allocate memory
      ReconstructionTypeMap = new char* [Ni];
      for (int i = 0; i < Ni ; ++i){
	ReconstructionTypeMap[i] = new char [Nj];
      }// endfor
    }
  } else {
    deallocate_ReconstructionTypeMap();
  }
  
}

// deallocate()
/*!
 * Deallocate all allocated memory.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::deallocate(void){

  int i;

  if (_allocated_block){

    for ( i = 0; i < Ni ; ++i ) {
      delete [] TD[i]; TD[i] = NULL;  // deallocate TD
      delete [] SI[i]; SI[i] = NULL;  // deallocate SI
      delete [] LimitedCell[i]; LimitedCell[i] = NULL; // deallocate monotonicity flag
      delete [] PreviousLimitedCell[i]; PreviousLimitedCell[i] = NULL; // deallocate monotonicity flag copy 

      if (_allocated_psinv){
	delete [] CENO_LHS[i]; CENO_LHS[i] = NULL; // deallocate pseudo-inverse
	delete [] CENO_Geometric_Weights[i]; CENO_Geometric_Weights[i] = NULL; // deallocate geometric weights
      }

      if (ReconstructionTypeMap != NULL){
	delete [] ReconstructionTypeMap[i]; ReconstructionTypeMap[i] = NULL; // deallocate ReconstructionTypeMap
      }
    }

    delete [] TD; TD = NULL;
    delete [] SI; SI = NULL;
    delete [] LimitedCell; LimitedCell = NULL;
    delete [] PreviousLimitedCell; PreviousLimitedCell = NULL;

    if (_allocated_psinv){
      delete [] CENO_LHS; CENO_LHS = NULL;
      delete [] CENO_Geometric_Weights; CENO_Geometric_Weights = NULL;
    }

    if (ReconstructionTypeMap != NULL){
      delete [] ReconstructionTypeMap; ReconstructionTypeMap = NULL;
    }

    // Reset all indexes
    Ni = 0; Nj = 0; Ng = 0;
    ICl = 0; ICu = 0; JCl = 0; JCu = 0;
    OrderOfReconstruction = -1;
    Nghost_HO = 0;
    rings = 0;
    rings_SI = 0;
    StartI = 0; EndI = 0;
    StartJ = 0; EndJ = 0;

    // Separate the high-order object from the associated geometry
    Geom = NULL;
    // Reset grid observers
    ObserverInteriorCellGeometryState = 0;
    ObserverGhostCellGeometryState = 0;
    ObserverCornerGhostCellGeometryState = 0;
    _constrained_block_reconstruction = false;

    // Reset block boundaries
    ResetBlockBoundaries();

    // Confirm the deallocation
    _allocated_block = false;
    _allocated_cells = false;
    _allocated_psinv = false;
    _calculated_psinv = false;

    // Set other flags
    _freeze_limiter = false;

  }
}

// deallocate_CellMemory()
/*! 
 * Deallocate all memory at the cell level.
 * Automatic deallocation is already provided when the 
 * variables are deleted, so there is no need to call
 * this function in the 'deallocate' routine.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::deallocate_CellMemory(void){

  int ii, jj;

  if (_allocated_cells){     // check that memory has been allocated

    for (jj  = JCl-Nghost_HO ; jj <= JCu+Nghost_HO ; ++jj ) {
      for ( ii = ICl-Nghost_HO ; ii <= ICu+Nghost_HO ; ++ii ) {  
  
	// deallocate TD
	TD[ii][jj].free_memory();

	// deallocate LimitedCell
	LimitedCell[ii][jj].clear();

	// deallocate PreviousLimitedCell
	PreviousLimitedCell[ii][jj].clear();
	
	// deallocate SmoothnessIndicator
	SI[ii][jj].clear();

	// deallocate LHS matrix and CENO geometric weights
	if (_allocated_psinv){
	  CENO_LHS[ii][jj].newsize(0,0);
	  CENO_Geometric_Weights[ii][jj].clear();
	}

      } /* endfor */
    } /* endfor */

    // reset flags 
    _calculated_psinv = false;

    // Deallocate reconstruction helper variables
    DeltaCellCenters.clear();
    i_index.clear();
    j_index.clear();
    A.newsize(0,0);
    GeometricWeights.clear();
    All_Delta_U.newsize(0,0);
    Delta_U.newsize(0);
    X.newsize(0);

    // Deallocate constrained reconstruction helper variables
    i_index_ave.clear();
    j_index_ave.clear();
    Constraints_Loc.clear();
    Constraints_Normals.clear();
    Constraints_BCs.clear();
    Approx_Constraints_Loc.clear();
    Approx_Constraints_Normals.clear();
    Approx_Constraints_BCs.clear();
    A_Assembled.newsize(0,0);
    All_U_Assembled.newsize(0,0);
    X_Assembled.newsize(0,0);
    

    // Confirm the deallocation
    OrderOfReconstruction = -1;
    Nghost_HO = 0;
    rings = 0;
    rings_SI = 0;
    StartI = 0; EndI = 0;
    StartJ = 0; EndJ = 0;

    _allocated_cells = false;
    _allocated_psinv = false;
  }//endif
}

/*!
 * Deallocate memory allocated to the ReconstructionTypeMap variable.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::deallocate_ReconstructionTypeMap(void){

  if (ReconstructionTypeMap != NULL){
    for (int i = 0; i < Ni ; ++i ) {
      delete [] ReconstructionTypeMap[i]; ReconstructionTypeMap[i] = NULL; // deallocate ReconstructionTypeMap
    }

    delete [] ReconstructionTypeMap; ReconstructionTypeMap = NULL;
  }
}

/*!
 * Reset information related to block boundaries in the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetBlockBoundaries(void){

  // West boundaries
  WestBnd.Reset();
  S_WestBnd.Reset();
  N_WestBnd.Reset();

  // East boundaries
  EastBnd.Reset();
  S_EastBnd.Reset();
  N_EastBnd.Reset();

  // North boundaries
  NorthBnd.Reset();
  W_NorthBnd.Reset();
  E_NorthBnd.Reset();

  // South boundaries
  SouthBnd.Reset();
  W_SouthBnd.Reset();
  E_SouthBnd.Reset();
}

/*!
 * Initialize completely the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeVariable(int ReconstructionOrder, GeometryType & Block,
						 const bool &_pseudo_inverse_allocation_){

  // Initialize the basic functionality of the high-order object
  InitializeBasicVariable(ReconstructionOrder, Block, _pseudo_inverse_allocation_);

  // Compute the pseudo-inverse if required
  ComputeReconstructionPseudoInverse();
}

/*! 
 * Initialize the basic functionality of the high-order object 
 * (i.e. memory allocation, indexes etc.).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeBasicVariable(int ReconstructionOrder, GeometryType & Block,
						      const bool &_pseudo_inverse_allocation_){

  // Allocate (re-allocate) memory for the high-order object.
  allocate(Block.ICu-Block.ICl+1,
	   Block.JCu-Block.JCl+1,
	   Block.Nghost,
	   _pseudo_inverse_allocation_,
	   ReconstructionOrder);

  // Associate geometry
  SetGeometryPointer(Block);

  // Set properties of block boundaries (i.e. which boundaries are constrained/unconstrained and which are opaque/transparent)
  SetPropertiesHighOrderBlock();

  // Build the reconstruction type map
  BuildReconstructionTypeMap();

  // Compute the smoothness indicator adjustment coefficient
  /* It adjusts the value of the smoothness indicator based on the degree of the reconstruction */
  AdjustmentCoeff = (getStencilSize() - NumberOfTaylorDerivatives() )/( NumberOfTaylorDerivatives() - 1.0);
}

//! Reset the reconstruction order.
//  This function assumes that the object has already associated geometry.
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetReconstructionOrder(int ReconstructionOrder){

  // Use the already associated grid and pseudo-inverse allocation flag to set the new reconstruction order
  bool _pseudo_inverse_allocation_ = _allocated_psinv;

  InitializeVariable(ReconstructionOrder,
		     *Geom,
		     _pseudo_inverse_allocation_);
}

/*! 
 * Set the properties of high-order block boundaries,
 * based on the associated grid.
 * There are three possible boundary types:
 * 1. constrained and opaque
 * 2. unconstrained and opaque
 * 3. unconstrained and transparent
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetPropertiesHighOrderBlock(void){

  // Ensure that block splines require a valid scenario for constrained reconstruction
  CheckConsistencyOfGeometricSplineProperties();

  // Set properties of boundaries
  
  // === Reset all boundary properties (i.e. all boundaries are set to be unconstrained and transparent)
  ResetBlockBoundaries();
  // === Set constrained block to false
  _constrained_block_reconstruction = false;

  // === North boundary ===
  if (Geom->IsNorthBoundaryReconstructionConstrained()){
    // Set North boundary to be constrained
    NorthBnd.setConstrainedBoundary();

    // Set the North extensions of West and East boundaries to be opaque (see Rule #1)
    N_WestBnd.setOpaqueBoundary();
    N_EastBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }

  // === West boundary ===
  if (Geom->IsWestBoundaryReconstructionConstrained()){
    // Set West boundary to be constrained
    WestBnd.setConstrainedBoundary();
    
    // Set the West extensions of North and South boundaries to be unconstrained and opaque (see Rule #1)
    W_NorthBnd.setOpaqueBoundary();
    W_SouthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }

  // === South boundary ===
  if (Geom->IsSouthBoundaryReconstructionConstrained()){
    // Set South boundary to be constrained
    SouthBnd.setConstrainedBoundary();

    // Set the South extensions of West and East boundaries to be unconstrained and opaque (see Rule #1)
    S_WestBnd.setOpaqueBoundary();
    S_EastBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }

  // === East boundary ===
  if (Geom->IsEastBoundaryReconstructionConstrained()){
    // Set East boundary to be constrained
    EastBnd.setConstrainedBoundary();
    
    // Set the East extensions of North and South boundaries to be unconstrained and opaque (see Rule #1)
    E_NorthBnd.setOpaqueBoundary();
    E_SouthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }

  // === North-West corner ===
  if (Geom->IsNorthExtendWestBoundaryReconstructionConstrained()){
    // Set North extension of West boundary to be constrained
    N_WestBnd.setConstrainedBoundary();
    // Set West extension of North boundary to be opaque (i.e. avoid corner ghost cells)
    W_NorthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }
  if (Geom->IsWestExtendNorthBoundaryReconstructionConstrained()){
    // Set West extension of North boundary to be constrained
    W_NorthBnd.setConstrainedBoundary();
    // Set North extension of West boundary to be opaque (i.e. avoid corner ghost cells)
    N_WestBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }


  // === South-West corner ===
  if (Geom->IsWestExtendSouthBoundaryReconstructionConstrained()){
    // Set West extension of South boundary to be constrained
    W_SouthBnd.setConstrainedBoundary();
    // Set South extension of West boundary to be opaque (i.e. avoid corner ghost cells)
    S_WestBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }
  if (Geom->IsSouthExtendWestBoundaryReconstructionConstrained()){
    // Set South extension of West boundary to be constrained
    S_WestBnd.setConstrainedBoundary();
    // Set West extension of South boundary to be opaque (i.e. avoid corner ghost cells)
    W_SouthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }


  // === South-East corner ===
  if (Geom->IsSouthExtendEastBoundaryReconstructionConstrained()){
    // Set South extension of East boundary to be constrained
    S_EastBnd.setConstrainedBoundary();
    // Set East extension of South boundary to be opaque (i.e. avoid corner ghost cells)
    E_SouthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }
  if (Geom->IsEastExtendSouthBoundaryReconstructionConstrained()){
    // Set East extension of South boundary to be constrained
    E_SouthBnd.setConstrainedBoundary();
    // Set South extension of East boundary to be opaque (i.e. avoid corner ghost cells)
    S_EastBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }


  // === North-East corner ===
  if (Geom->IsEastExtendNorthBoundaryReconstructionConstrained()){
    // Set East extension of North boundary to be constrained
    E_NorthBnd.setConstrainedBoundary();
    // Set North extension of East boundary to be opaque (i.e. avoid corner ghost cells)
    N_EastBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }
  if (Geom->IsNorthExtendEastBoundaryReconstructionConstrained()){
    // Set North extension of East boundary to be constrained
    N_EastBnd.setConstrainedBoundary();
    // Set East extension of North boundary to be opaque (i.e. avoid corner ghost cells)
    E_NorthBnd.setOpaqueBoundary();

    // Flag block as constrained
    _constrained_block_reconstruction = true;
  }

  // Set indexes between which solution reconstruction is performed.
  StartI = ICl - Nghost_HO; EndI = ICu + Nghost_HO;
  StartJ = JCl - Nghost_HO; EndJ = JCu + Nghost_HO;

}

/*! 
 * Ensure that the block geometric splines require a realizable
 * constrained reconstruction scenario (i.e. some rules are fulfilled).
 * 
 * Rule 1. If one main spline is constrained, then the extensions of the main splines
 *         perpendincular on it and which have contact with it, cannot be constrained.
 *         (e.g. North spline constrained, then ExtendNorth_WestSpline and 
 *               ExtendNorth_EastSpline cannot be constrained.)
 * Rule 2. If both main splines that meet in a corner are constrained,
 *         then the extension splines at that corner cannot be constrained. 
 * Rule 3. If one extension spline is constrained, the other extension spline in the corner
 *         and its correspondent main spline cannot be constrained 
 *         (e.g. ExtendEast_SouthSpline constrained, then ExtendSouth_EastSpline and EastSpline
 *          cannot be constrained)
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties(void){

  // === Check main North boundary ===
  if (Geom->IsNorthBoundaryReconstructionConstrained()){
    // === Check Rule 1 ===
    if (Geom->IsNorthExtendWestBoundaryReconstructionConstrained() ||
	Geom->IsNorthExtendEastBoundaryReconstructionConstrained() ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #1 violated for the constrained North boundary");
    }

    // === Check Rule 2 for NW corner ===
    if (Geom->IsWestBoundaryReconstructionConstrained() && 
	(Geom->IsWestExtendNorthBoundaryReconstructionConstrained() || Geom->IsNorthExtendWestBoundaryReconstructionConstrained()) ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #2 violated for the North-West corner");
    }

    // === Check Rule 2 for NE corner ===
    if (Geom->IsEastBoundaryReconstructionConstrained() && 
	(Geom->IsEastExtendNorthBoundaryReconstructionConstrained() || Geom->IsNorthExtendEastBoundaryReconstructionConstrained()) ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #2 violated for the North-East corner");
    }
  }

  // === Check main South boundary ===
  if (Geom->IsSouthBoundaryReconstructionConstrained()){
    // === Check Rule 1 ===
    if (Geom->IsSouthExtendWestBoundaryReconstructionConstrained() ||
	Geom->IsSouthExtendEastBoundaryReconstructionConstrained() ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #1 violated for the constrained South boundary");
    }

    // === Check Rule 2 for SW corner ===
    if (Geom->IsWestBoundaryReconstructionConstrained() && 
	(Geom->IsWestExtendSouthBoundaryReconstructionConstrained() || Geom->IsSouthExtendWestBoundaryReconstructionConstrained()) ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #2 violated for the South-West corner");
    }

    // === Check Rule 2 for SE corner ===
    if (Geom->IsEastBoundaryReconstructionConstrained() && 
	(Geom->IsEastExtendSouthBoundaryReconstructionConstrained() || Geom->IsSouthExtendEastBoundaryReconstructionConstrained()) ){
      throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #2 violated for the South-East corner");
    }
  }

  // === Check Rule 1 for main East boundary (Obs: Rule #2 has been already checked) ===
  if (Geom->IsEastBoundaryReconstructionConstrained() && (Geom->IsEastExtendNorthBoundaryReconstructionConstrained() || 
							  Geom->IsEastExtendSouthBoundaryReconstructionConstrained()) ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #1 violated for the constrained East boundary");
  }

  // === Check Rule 1 for main West boundary (Obs: Rule #2 has been already checked) ===
  if (Geom->IsWestBoundaryReconstructionConstrained() && (Geom->IsWestExtendNorthBoundaryReconstructionConstrained() || 
							  Geom->IsWestExtendSouthBoundaryReconstructionConstrained()) ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #1 violated for the constrained West boundary");
  }

  // === Check Rule 3 for North-West corner (Obs: It's enough to check only the two extensions) ===
  if (Geom->IsWestExtendNorthBoundaryReconstructionConstrained() && 
      Geom->IsNorthExtendWestBoundaryReconstructionConstrained() ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #3 violated for the North-West corner");
  }

  // === Check Rule 3 for North-East corner (Obs: It's enough to check only the two extensions) ===
  if (Geom->IsEastExtendNorthBoundaryReconstructionConstrained() && 
      Geom->IsNorthExtendEastBoundaryReconstructionConstrained() ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #3 violated for the North-East corner");
  }

  // === Check Rule 3 for South-West corner (Obs: It's enough to check only the two extensions) ===
  if (Geom->IsWestExtendSouthBoundaryReconstructionConstrained() && 
      Geom->IsSouthExtendWestBoundaryReconstructionConstrained() ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #3 violated for the South-West corner");
  }

  // === Check Rule 3 for South-East corner (Obs: It's enough to check only the two extensions) ===
  if (Geom->IsEastExtendSouthBoundaryReconstructionConstrained() && 
      Geom->IsSouthExtendEastBoundaryReconstructionConstrained() ){
    throw runtime_error("HighOrder2D<SOLN_STATE>::CheckConsistencyOfGeometricSplineProperties() ERROR! Consistency rule #3 violated for the South-East corner");
  }
}

// AssociateGeometry()
/*! 
 * Perform all settings when a specific geometry is
 * associated to the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::AssociateGeometry(GeometryType & Block){

  bool _pseudo_inverse_allocation_ = _allocated_psinv;

  InitializeVariable(OrderOfReconstruction,
		     Block,
		     _pseudo_inverse_allocation_);
}

// getNumberOfRings()
/*! Return the required number of neighbor rings
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::getNumberOfRings(int number_of_Taylor_derivatives){
  switch(number_of_Taylor_derivatives){
  case 1:   // piecewise constant
    return 0;			// it corresponds to 0 neighbour cells

  case 3:   // piecewise linear
    return 1; 			// it corresponds to 8 neighbour cells

  case 6:   // piecewise quadratic
    return 2;			// it corresponds to 24 neighbour cells

  case 10:  // piecewise cubic
    return 2;			// it corresponds to 24 neighbour cells

  case 15:  // piecewise quartic
    return 2;			// it corresponds to 24 neighbour cells

  case 21:  // piecewise quintic
    return 3;			// it corresponds to 48 neighbour cells

  default:
    return 0;
  }
}

// getTaylorDerivativesSize()
/*! Get the number of Taylor derivatives based on
 * the reconstruction order.
 * In 2D, this number is given by \f$ \frac{(k+1) \, (k+2)}{2} \f$,
 * where k is the order of reconstruction.
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::getTaylorDerivativesSize(const int &ReconstructionOrder) {
  return (ReconstructionOrder + 1) * (ReconstructionOrder + 2)/2;
}

// getNghostHighOrder()
/*! 
 * Return the number of ghost cells in which high-order reconstruction
 * needs to be performed based on the provided reconstruction order.
 * This number takes into account factors such as: \n
 * - number of ghost cells required for flux calculation (i.e. 1). \n
 * - stencil size for a given order of reconstruction. \n
 * - number of cells required to compute the smoothness indicator. \n
 * - flags set in the Execution_Mode class.
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::getNghostHighOrder(const int &ReconstructionOrder){

  if (CENO_Execution_Mode::CENO_RECONSTRUCTION_WITH_MESSAGE_PASSING){
    // This method is not implemented yet but it's left as a possible development stream.
    // In the current implementation the 'CENO_Execution_Mode' class should refuse to set this flag to ON.

  } else {

    if (CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS){
      // Return the minimum number of ghost cells when smoothness indicator computation is done with only the closes neighbours
      switch (ReconstructionOrder){
      case 4:
	return 2;
	break;

      case 3:
	return 2;
	break;

      case 2:
	return 2;
	break;

      case 1:
	return 2;
	break;

      case 0:
	return 1;
	break;

      case -1:
	return 0;
	break;

      default:
	// require a very large number such that to generate a stop
	return 20;
      }	// endswitch

    } else {
      // Return the minimum number of ghost cells when smoothness indicator computation is done with all used neighbours
      switch (ReconstructionOrder){
      case 4:
	return 3;
	break;

      case 3:
	return 3;
	break;

      case 2:
	return 3;
	break;

      case 1:
	return 2;
	break;

      case 0:
	return 1;
	break;

      case -1:
	return 0;
	break;

      default:
	// require a very large number such that to generate a stop
	return 20;
      }	// endswitch
    }

  }//endif    
}

// getMinimumNghost()
/*!
 * Returns the minimum number of ghost cells for the solution domain
 * required to carry out a reconstruction of order ReconstructionOrder.
 * This number takes into account factors such as: \n
 * - number of ghost cells required for flux calculation (i.e. 1). \n
 * - stencil size for a given order of reconstruction. \n
 * - number of cells required to compute the smoothness indicator. \n
 * - flags set in the Execution_Mode class.
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::getMinimumNghost(const int &ReconstructionOrder){
  return getNghostHighOrder(ReconstructionOrder) + getNumberOfRings(getTaylorDerivativesSize(ReconstructionOrder));
}

// SetRings()
/*! Set the number of rings around the current cell
 * which will be used to form the supporting stencil
 * of the reconstruction. 
 * The number of derivatives is determined with 
 * getTaylorDerivativesSize() subroutine based on
 * the order of reconstruction.
 * Set also the number of rings around the current cell
 * which will be used to compute the smoothness indicator.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetRings(void){
  rings = getNumberOfRings(getTaylorDerivativesSize());

  if (CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS){
    // use only one layer of neighbours
    rings_SI = 1;
  } else {
    // use the same number of rings as the reconstruction
    rings_SI = rings;
  }
}

//!< getStencilSize()
/*!
 * Return the stencil size for a CENO reconstruction of a given reconstruction order.
 */
template<class SOLN_STATE> inline
int HighOrder2D<SOLN_STATE>::getStencilSize(const int &ReconstructionOrder) {
  static int temp;

  // Calculate temp based on the number of rings
  temp = 1 + 2*getNumberOfRings(getTaylorDerivativesSize(ReconstructionOrder));

  return temp*temp;
} 

//! Reset the monotonicity data (i.e. flag + limiter)
//  throughout the block.
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetMonotonicityData(void){

  int i,j;

  for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
    for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i ) {    
      ResetMonotonicityData(i,j);
    } /* endfor */
  } /* endfor */

}

/*!
 * Reset the monotonicity data (i.e. flag + limiter)
 * for the specified cell.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetMonotonicityData(const int & ii, const int & jj){

  /* Reset the CENO smoothness indicator flag to OFF (i.e. smooth solution) */
  for (int k = 0; k <= NumberOfVariables() - 1; ++k){
    LimitedCell[ii][jj][k] = OFF;
  } /* endfor */
  
  // Reset the limiter value (i.e. no limiting)
  TD[ii][jj].ResetLimiter(); 
}

/*! Reset the backup of monotonicity data (i.e. flag + limiter)
 *  throughout the block. 
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetMonotonicityDataBackup(void){

  int i,j;

  for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
    for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i ) {    
      ResetMonotonicityDataBackup(i,j);
    } /* endfor */
  } /* endfor */

}

/*!
 * Reset the backup of monotonicity data 
 * (i.e. flag copy + frozen limiter) for the specified cell.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetMonotonicityDataBackup(const int & ii, const int & jj){

  /* Reset the CENO smoothness indicator flag copy to OFF (i.e. smooth solution) */
  for (int k = 0; k <= NumberOfVariables() - 1; ++k){
    PreviousLimitedCell[ii][jj][k] = OFF;
  } /* endfor */
  
  // Reset the value of the frozen limiter (i.e. no limiting)
  TD[ii][jj].ResetFrozenLimiter();
}

/*! 
 * Allocate memory and initialize the variables 
 * used for storing monotonicity information.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeMonotonicityVariables(const int & ii, const int & jj){

  // Allocate memory and initialize containers
  SI[ii][jj].assign(NumberOfVariables(), 0.0);
  LimitedCell[ii][jj].assign(NumberOfVariables(), OFF);
  PreviousLimitedCell[ii][jj].assign(NumberOfVariables(), OFF);

}

/*! 
 * Assess the smoothness of each solution variable interpolant
 * relative to the CENO_Tolerances::Fit_Tolerance.
 * Flag with inadequate reconstructions those interpolants detected
 * as non-smooth.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::AssessInterpolantsSmoothness(const int &iCell, const int &jCell){

  int parameter;

  for(parameter=1; parameter<=NumberOfVariables(); ++parameter){

    if( CellSmoothnessIndicatorValue(iCell,jCell,parameter) < CENO_Tolerances::Fit_Tolerance ){

      // Flag the (iCell,jCell) cell with inadequate reconstruction for the current variable
      CellInadequateFitValue(iCell,jCell,parameter) = ON;

      if (CENO_Execution_Mode::CENO_PADDING){
	/* Flag also all cells surrounding the (iCell,jCell) cell with inadequate reconstruction if CENO_Padding is ON */
	CellInadequateFitValue(iCell-1,jCell-1,parameter) = ON;
	CellInadequateFitValue(iCell  ,jCell-1,parameter) = ON;
	CellInadequateFitValue(iCell+1,jCell-1,parameter) = ON;

	CellInadequateFitValue(iCell-1,jCell ,parameter) = ON;
	CellInadequateFitValue(iCell+1,jCell ,parameter) = ON;

	CellInadequateFitValue(iCell-1,jCell+1,parameter) = ON;
	CellInadequateFitValue(iCell  ,jCell+1,parameter) = ON;
	CellInadequateFitValue(iCell+1,jCell+1,parameter) = ON;
      }//endif
    }//endif

  }//endfor(parameter)

}

/*! 
 * Compute the CENO smoothness indicator which is used to
 * differentiate between smooth and non-smooth solution content
 * for all block cells which are used for flux calculation.
 * Flag also those cells detected with non-smooth solution content.
 * The typical range of cells if formed by the interior cells and
 * the first layer of ghost cells. 
 * If constrained reconstruction is used along some boundaries the
 * ghost cells adjacent to those don't need to have the indicator computed.
 *
 * \param [in] SolnBlk The solution block which provides solution data.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the average solution.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk,
							 const Soln_State &
							 (Soln_Block_Type::*ReconstructedSoln)(const int &,
											       const int &) const){

  // Set local variables
  int parameter;
  int i,j;
  
  // == Check if the reconstruction polynomial is piecewise constant
  if (RecOrder() == 0){
    // There is no need to calculate the smoothness indicator
    return;
  }

  
  // Check whether constrained reconstruction is required anywhere in the block.
  if ( !_constrained_block_reconstruction ){

    // ===  Compute the smoothness indicator using the central stencil ===
    
    // Compute the smoothness indicator for cells inside the domain
    for ( j  = JCl ; j <= JCu ; ++j ) {
      for ( i = ICl ; i <= ICu ; ++i ) {
	ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						    ReconstructedSoln,
						    i, j);	
      } /* endfor(i) */
    } /* endfor(j) */
    
    // Compute the smoothness indicator for ghost cells involved in flux calculation
    // == South and North boundaries ==
    for (i = ICl; i<=ICu; ++i){
      // == South bnd.
      ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						  ReconstructedSoln,
						  i, JCl-1);	
      
      // == North bnd.
      ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						  ReconstructedSoln,
						  i, JCu+1);
    }

    // == West and East boundaries ==
    for (j = JCl; j<=JCu; ++j){
      // == West bnd.
      ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						  ReconstructedSoln,
						  ICl-1, j);      

      // == East bnd.
      ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						  ReconstructedSoln,
						  ICu+1, j);      
    }

  } else {

    // ===  Compute the smoothness indicator using the central and deviated stencils ===

    // Compute the smoothness indicator for cells inside the domain
    for ( j  = JCl ; j <= JCu ; ++j ) {
      for ( i = ICl ; i <= ICu ; ++i ) {
	ComputeCellSmoothnessDataForConstrainedBlock(SolnBlk,
						     ReconstructedSoln,
						     i,j);
      } /* endfor(i) */
    } /* endfor(j) */    

    // Compute the smoothness indicator for ghost cells involved in flux calculation
    // === West Bnd. ===
    if (!IsWestConstrainedReconstructionRequired()){
      for (j = JCl; j <= JCu; ++j){
	ComputeCellSmoothnessDataForConstrainedBlock(SolnBlk,
						     ReconstructedSoln,
						     ICl-1,j);
      }
    }

    // === East Bnd. ===
    if (!IsEastConstrainedReconstructionRequired()){
      for (j = JCl; j <= JCu; ++j){
	ComputeCellSmoothnessDataForConstrainedBlock(SolnBlk,
						     ReconstructedSoln,
						     ICu+1,j);
      }
    }

    // === North Bnd. ===
    if (!IsEastConstrainedReconstructionRequired()){
      for (i = ICl; i <= ICu; ++i){
	ComputeCellSmoothnessDataForConstrainedBlock(SolnBlk,
						     ReconstructedSoln,
						     i,JCu+1);
      }
    }

    // === South Bnd. ===
    if (!IsEastConstrainedReconstructionRequired()){
      for (i = ICl; i <= ICu; ++i){
	ComputeCellSmoothnessDataForConstrainedBlock(SolnBlk,
						     ReconstructedSoln,
						     i,JCl-1);
      }
    }

  } // endif (_constrained_block_reconstruction)

}

/*! 
 * Compute the data related to CENO smoothness indicator
 * for a given cell (iCell,jCell) which is know to have a 
 * central stencil.
 *
 * \param [in] SolnBlk The solution block which provides solution data.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the average solution.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeCellSmoothnessDataWithCentralStencil(Soln_Block_Type &SolnBlk,
									  const Soln_State & 
									  (Soln_Block_Type::*ReconstructedSoln)(const int &,
														const int &) const,
									  const int &iCell, const int &jCell){


  /* 
   * Set the central supporting stencil.
   * Write the 'i' and 'j' indexes of the cells that are part of
   * the smoothness indicator calculation stencil of cell (iCell,jCell).
   * Use the number of smoothness indicator rings set in the class to 
   * determine how far the stencil extends.
   * The indexes are written in the i_index and j_index containers.
   * Because these containers might be larger than needed, 
   * the StencilSize_SmoothnessIndicator variable is used for storing 
   * the smoothness indicator stencil size.
   *
   * Note: The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
   */
  SetCentralStencil(iCell,jCell,
		    i_index,j_index,
		    rings_SI,StencilSize_SmoothnessIndicator);
	
  // Evaluate the Smoothness Indicator for the current cell for all solution state variables
  ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
			     iCell, jCell, i_index, j_index, StencilSize_SmoothnessIndicator);
  
  // Check the smoothness condition
  AssessInterpolantsSmoothness(iCell,jCell);

}

/*! 
 * Compute the data related to CENO smoothness indicator
 * for a given cell (iCell,jCell) which is know to have a 
 * deviated stencil.
 *
 * \param [in] SolnBlk The solution block which provides solution data.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the average solution.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeCellSmoothnessDataWithDeviatedStencil(Soln_Block_Type &SolnBlk,
									   const Soln_State &
									   (Soln_Block_Type::*ReconstructedSoln)(const int &,
														 const int &) const,
									   const int &iCell, const int &jCell){
  
  // Set the biased supporting stencil (i.e. deviated from central and not extended)
  SetDeviatedReconstructionStencil(iCell,jCell,
				   i_index_ave, j_index_ave,
				   RingsSI(), false);

  // Evaluate the Smoothness Indicator for the current cell for all solution state variables
  ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
			     iCell, jCell, i_index_ave, j_index_ave, i_index_ave.size());
  
  // Check the smoothness condition
  AssessInterpolantsSmoothness(iCell,jCell);

}

/*! 
 * Compute the data related to CENO smoothness indicator
 * for a given cell (iCell,jCell) which is know to be 
 * in a constrained block.
 *
 * The purpose of this function is to provide a compact routine which 
 * picks the proper method for computing smoothness indicator data
 * based on the reconstruction type of the given cell.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeCellSmoothnessDataForConstrainedBlock(Soln_Block_Type &SolnBlk,
									   const Soln_State &
									   (Soln_Block_Type::*ReconstructedSoln)(const int &,
														 const int &) const,
									   const int &iCell, const int &jCell){
  
  switch (ReconstructionTypeMap[iCell][jCell]){
  case 'r':		// "Regular reconstruction"
    // Compute the smoothness indicator using the central stencil
    ComputeCellSmoothnessDataWithCentralStencil(SolnBlk,
						ReconstructedSoln,
						iCell, jCell);
    break;

  case 'm':		// "Modified reconstruction" 
  case 'c':		// "Constrained reconstruction"
    // Compute the smoothness indicator using the deviated stencil
    ComputeCellSmoothnessDataWithDeviatedStencil(SolnBlk,
						 ReconstructedSoln,
						 iCell, jCell);
    break;
	  
  case 'n':		// "No reconstruction"
    // Do nothing
    break;
  } //endswitch

}

/*! 
 * Compute the CENO smoothness indicator which is used to
 * differentiate between smooth and non-smooth solution content
 * for a specific cell (iCell,jCell).
 *
 * \param [in] SolnBlk The solution block which provides solution data.
 * \param ReconstructedSoln member function of Soln_Block_Type which returns the average solution.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk,
							 const Soln_State &
							 (Soln_Block_Type::*ReconstructedSoln)(const int &,
											       const int &) const,
							 const int &iCell, const int &jCell,
							 const IndexType & i_index, const IndexType & j_index,
							 const int & StencilSize){

  // SET VARIABLES USED IN THE ANALYSIS PROCESS
  double alpha;
  int parameter, cell;

  // Set the mean solution of (iCell,jCell). It's used as a reference.
  MeanSolution = (SolnBlk.*ReconstructedSoln)(iCell,jCell);

  // Initialize the minimum mean solution in absolute value.
  Min_Mean_Solution = fabs(MeanSolution);

  /*! NOTE: The following used variables are set as private to the class:
    SS_Regression, SS_Residual, Max_SS_Regression, Temp_Regression, Temp_Residual. */

  // Initialize SS_Regression and SS_Residual
  SS_Regression  = CellTaylorDerivState(iCell,jCell,0) - MeanSolution;
  SS_Regression *= SS_Regression;      // get the squared value
  Max_SS_Regression = SS_Regression;   // initialize Max_SS_Regression 
  SS_Residual.Vacuum(); // Note: the residual difference for (iCell,jCell) is zero.

  // Get the normalization state
  NormalizationState = SolnBlk.getNormalizationState(iCell,jCell);

  /* Compute SS_Regression and SS_Residual for all the variables at once */
  for( cell=1; cell<StencilSize; ++cell){

    // Get minimum average solution in absolute value
    Min_Mean_Solution = min(Min_Mean_Solution , fabs((SolnBlk.*ReconstructedSoln)(i_index[cell],j_index[cell])));

    Temp_Regression  = CellTaylorDerivState(i_index[cell],j_index[cell],0) - MeanSolution;
    Temp_Regression *= Temp_Regression;          /* compute Temp_Regression square */

    // Get maximum squared solution variation for the current Temp_Regression
    Max_SS_Regression = max(Max_SS_Regression, Temp_Regression);

    Temp_Residual  = ( CellTaylorDerivState(i_index[cell],j_index[cell],0) - 
		       SolutionStateAtLocation(iCell,jCell,CellCenter(i_index[cell],j_index[cell])) );

    // Update SS_Regression & SS_Residual
    if (CENO_Execution_Mode::CENO_CONSIDER_WEIGHTS){

      /* Compute the X and Y component of the distance between
	 the cell center of the neighbour cell and the reconstructed cell */
      DeltaCentroids = CellCenter(i_index[cell],j_index[cell]) - CellCenter(iCell,jCell);
    
      // Compute the geometric weight based on the centroid distance
      CENO_Geometric_Weighting(GeomWeightSI, DeltaCentroids.abs());

      // Update SS_Regression
      SS_Regression += Temp_Regression * GeomWeightSI;
      // Update SS_Residual
      Temp_Residual *= Temp_Residual;       /* compute Temp_Residual square */
      SS_Residual   += Temp_Residual * GeomWeightSI;

    } else {

      // Update SS_Regression
      SS_Regression += Temp_Regression;

      // Update SS_Residual
      Temp_Residual *= Temp_Residual;       /* compute Temp_Residual square */
      SS_Residual  += Temp_Residual;
    }//endif

  }//endfor(cell)

  /* Compute the smoothness indicator for each variable */
  for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){

    /* Decide if the 'alpha' for the current parameter is computed or not.
       This decision is dictated by the following reasons:
       --> If there is not at all or very small solution variation within the stencil (i.e. uniform flow)
           SS_Residual[parameter] is approximately equal to SS_Regression[parameter].
	   That would trigger 'alpha' to have a very small value and consequently the cell flagged as unfit.
       --> The user should have the freedom to specify a level of numerical noise under which solution
           discontinuities are not of interest.
       To solve both these problems, the maximum squared variation is compared relative to an accepted 
       level of solution variation.
       To obtain consistency in setting the tolerated lack of smoothness for large ranges of solution values,
       a normalization is employed with the reference state provided by the solution state class for the 
       current parameter.
    */
    if ( Max_SS_Regression[parameter]/sqr(NormalizationState[parameter]) > 
	 CENO_Tolerances::SquareToleranceAroundValue(Min_Mean_Solution[parameter]/NormalizationState[parameter]) ){
      
      // Compute 'alpha'
      alpha = 1.0 - (SS_Residual[parameter]/SS_Regression[parameter]);
      
    } else {
      
      // There is not enough variation in the solution based on user set tolerance to flag a potential discontinuity.
      // Assign the perfect fit value to the smoothness indicator
      alpha = 1.0;
      
    } // endif

    /* Compute final value */
    CellSmoothnessIndicatorValue(iCell,jCell,parameter) = (alpha/(max(CENO_Tolerances::epsilon,1.0 - alpha)))*AdjustmentCoeff;

  } // endfor (parameter)

}

/*! 
 * Return true if any of the high-order solution reconstructions
 * of cell (iCell,jCell) is flagged as non-smooth, otherwise
 * return false.
 */
template<class SOLN_STATE> inline
const bool HighOrder2D<SOLN_STATE>::IsThereAnyNonSmoothHighOrderReconstruction(const int &iCell,
									       const int &jCell) const {

  // Analyse 'CellInadequateFit' flag for each solution variable.
  for(int parameter = 1; parameter <= NumberOfVariables(); ++parameter){
    if ( CellInadequateFitValue(iCell,jCell,parameter) ){
      return true; 	// found already one interpolant flagged as non-smooth
    }//endif
  }//endfor(parameter)

  // No interpolant flagged as non-smooth was found.
  return false;
}

/*! 
 * This function bypass the smoothness indicator analysis
 * and flags all reconstructions of (iCell,jCell) cell 
 * as non-smooth.
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::FlagCellReconstructionsAsNonSmooth(const int &iCell, const int &jCell){
  for(int parameter = 1; parameter <= NumberOfVariables(); ++parameter){
    CellInadequateFitValue(iCell,jCell,parameter) = ON;
  }//endfor(parameter)  
}

/*!
 * Evaluate the gradient in the normal direction.
 * \param norm_dir Cartesian vector that defines the direction in which the gradient is computed.
 */
template<class SOLN_STATE> inline
SOLN_STATE HighOrder2D<SOLN_STATE>::NormalGradientStateAtCoordinates(const int & ii, const int & jj,
								     const double &X_Coord, const double &Y_Coord,
								     const Vector2D &norm_dir) const {
  return ( norm_dir.x*XGradientStateAtCoordinates(ii,jj,X_Coord,Y_Coord) + 
	   norm_dir.y*YGradientStateAtCoordinates(ii,jj,X_Coord,Y_Coord) );
}

/*!
 * Evaluate the gradient in the normal direction for a specific variable.
 * \param norm_dir Cartesian vector that defines the direction in which the gradient is computed.
 * \param parameter variable index.
 */
template<class SOLN_STATE> inline
double HighOrder2D<SOLN_STATE>::NormalGradientAtCoordinates(const int & ii, const int & jj,
							    const double &X_Coord, const double &Y_Coord,
							    const Vector2D &norm_dir, const unsigned & parameter) const {
  
  return ( norm_dir.x * TD[ii][jj].ComputeXGradientFor(X_Coord - XCellCenter(ii,jj),
						       Y_Coord - YCellCenter(ii,jj),
						       parameter) +
	   norm_dir.y * TD[ii][jj].ComputeYGradientFor(X_Coord - XCellCenter(ii,jj),
						       Y_Coord - YCellCenter(ii,jj),
						       parameter) );
}

//! Integrate over the cell
template<class SOLN_STATE>
template<typename FO, class ReturnType> inline
ReturnType HighOrder2D<SOLN_STATE>::IntegrateOverTheCell(const int &ii, const int &jj,
							 const FO FuncObj,
							 const int & digits,
							 ReturnType _dummy_param) const {
  return Geom->Integration.IntegrateFunctionOverCell(ii,jj, FuncObj, digits, _dummy_param);
}

/*! 
 * Integrate the reconstruction of cell (ii,jj)
 * over the specified quadrilateral domain.
 *
 * \param ii i-index of the cell
 * \param jj j-index of the cell
 * \param SW the South-West quadrilateral vertex
 * \param SE the South-East quadrilateral vertex
 * \param NE the North-East quadrilateral vertex
 * \param NW the North-West quadrilateral vertex
 */
template<class SOLN_STATE>
template<typename Node2DType> inline
SOLN_STATE HighOrder2D<SOLN_STATE>::
IntegrateCellReconstructionOverQuadDomain(const int &ii, const int &jj,
					  const Node2DType &SW, const Node2DType &NW,
					  const Node2DType &NE, const Node2DType &SE) const {

  SOLN_STATE _dummy_result;

  return ( Geom->Integration.
	   IntegratePolynomialOverQuadDomain(wrapped_soln_block_member_function(this,
										&ClassType::SolutionStateAtCoordinates,
										ii, jj,
										_dummy_result),
					     SW, NW, NE, SE,
					     _dummy_result) );
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the CENTRAL stencil of cell (iCell,jCell). To decide how far 
 * this stencil extends the routine uses the passed number of rings.
 *
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 * \param [in]  rings The number of neighbour cell rings around (iCell,jCell) cell.
 * \param [out] StencilSize This variable gets set to the stencil size.
 * This parameter is different than the class variable "rings".
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetCentralStencil(const int &iCell, const int &jCell,
						IndexType & i_index, IndexType & j_index,
						const int &rings, int &StencilSize) const{

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    i_index[9] =iCell-2;  j_index[9]=jCell-2;
    i_index[10]=iCell-1; j_index[10]=jCell-2;
    i_index[11]=iCell  ; j_index[11]=jCell-2;
    i_index[12]=iCell+1; j_index[12]=jCell-2;
    i_index[13]=iCell+2; j_index[13]=jCell-2;
    i_index[14]=iCell-2; j_index[14]=jCell-1;
    i_index[15]=iCell+2; j_index[15]=jCell-1;
    i_index[16]=iCell-2; j_index[16]=jCell;
    i_index[17]=iCell+2; j_index[17]=jCell;
    i_index[18]=iCell-2; j_index[18]=jCell+1;
    i_index[19]=iCell+2; j_index[19]=jCell+1;
    i_index[20]=iCell-2; j_index[20]=jCell+2;
    i_index[21]=iCell-1; j_index[21]=jCell+2;
    i_index[22]=iCell  ; j_index[22]=jCell+2;
    i_index[23]=iCell+1; j_index[23]=jCell+2;
    i_index[24]=iCell+2; j_index[24]=jCell+2;

  case 1: // one ring of cells around (iCell,jCell)

    i_index[0]=iCell;   j_index[0]=jCell; /* cell (iCell,jCell) */
    /* First ring */
    i_index[1]=iCell-1; j_index[1]=jCell-1;
    i_index[2]=iCell;   j_index[2]=jCell-1;
    i_index[3]=iCell+1; j_index[3]=jCell-1;
    i_index[4]=iCell-1; j_index[4]=jCell;
    i_index[5]=iCell+1; j_index[5]=jCell;
    i_index[6]=iCell-1; j_index[6]=jCell+1;
    i_index[7]=iCell;   j_index[7]=jCell+1;
    i_index[8]=iCell+1; j_index[8]=jCell+1;

    if (rings == 2){
      StencilSize = 25;
    } else {
      StencilSize = 9;
    }
    break;

  default: // general expression
    i_index[0] = iCell;
    j_index[0] = jCell;
    for (int i=iCell-rings, Poz=1; i<=iCell+rings; ++i){
      for (int j=jCell-rings; j<=jCell+rings; ++j){
	if(!((i==iCell)&&(j==jCell)) ){
	  i_index[Poz] = i;
	  j_index[Poz] = j;
	  ++Poz;
	}
      }
      StencilSize = Poz;
    }
  }//endswitch
  
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the reconstruction of cell (iCell,jCell). Use the number of
 * rings set in the class to determine how far the stencil extends.
 * This routine doesn't modify the stencil due to existence 
 * of curved boundaries.
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetReconstructionStencil(const int &iCell, const int &jCell,
						       IndexType & i_index, IndexType & j_index) const{

  int _dummy_;

  // Call set central stencil
  SetCentralStencil(iCell,jCell,i_index,j_index,rings,_dummy_);
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the reconstruction of cell (iCell,jCell) for special tests.
 * Use the number of rings set in the class to determine how far the stencil extends.
 * This routine doesn't modify the stencil due to existence 
 * of curved boundaries.
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetSpecialReconstructionStencil(const int &iCell, const int &jCell,
							      IndexType & i_index, IndexType & j_index) const{
  
  i_index.clear();
  j_index.clear();

  i_index.push_back(iCell);   j_index.push_back(jCell); /* cell (iCell,jCell) */

  switch(rings){

  case 2: // two rings of cells around (iCell,jCell)

    /* Second ring */
    //    i_index.push_back(iCell-2); j_index.push_back(jCell-2);
    i_index.push_back(iCell-1); j_index.push_back(jCell-2);
    //    i_index.push_back(iCell  ); j_index.push_back(jCell-2);
    i_index.push_back(iCell+1); j_index.push_back(jCell-2);
    //    i_index.push_back(iCell+2); j_index.push_back(jCell-2);
    i_index.push_back(iCell-2); j_index.push_back(jCell-1);
    i_index.push_back(iCell+2); j_index.push_back(jCell-1);
    //    i_index.push_back(iCell-2); j_index.push_back(jCell  );
    //    i_index.push_back(iCell+2); j_index.push_back(jCell  );
    i_index.push_back(iCell-2); j_index.push_back(jCell+1);
    i_index.push_back(iCell+2); j_index.push_back(jCell+1);
    //    i_index.push_back(iCell-2); j_index.push_back(jCell+2);
    i_index.push_back(iCell-1); j_index.push_back(jCell+2);
    //    i_index.push_back(iCell  ); j_index.push_back(jCell+2);
    i_index.push_back(iCell+1); j_index.push_back(jCell+2);
    //    i_index.push_back(iCell+2); j_index.push_back(jCell+2);

  case 1: // one ring of cells around (iCell,jCell)

    /* First ring */
    i_index.push_back(iCell-1); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell-1);
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell+1); j_index.push_back(jCell+1);

    break;

  default: // general expression

    break;
  }//endswitch  
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the reconstruction of the given cell and the cells that have
 * common face with it.
 *
 * \param [in]  iCell The i-index of the given cell
 * \param [in]  jCell The j-index of the given cell
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 *
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::getEnlargedReconstructionStencil(const int &iCell, const int &jCell,
							       IndexType & i_index, IndexType & j_index) const{

  if ( IsConstrainedReconstructionRequired() ) {
    throw runtime_error("HighOrder2D<SOLN_STATE>::getEnlargedReconstructionStencil() doesn't support constrained reconstruction!");
  }

  // Set stencil based on the central one

  // Reset the output index array
  i_index.clear();
  j_index.clear();

  switch(rings){

  case 2: 

    i_index.push_back(iCell);   j_index.push_back(jCell); /* cell (iCell,jCell) */

    /* First ring */
    i_index.push_back(iCell-1); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell-1);
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell+1); j_index.push_back(jCell+1);

    /* Second ring */
    i_index.push_back(iCell-2); j_index.push_back(jCell-2);
    i_index.push_back(iCell-1); j_index.push_back(jCell-2);
    i_index.push_back(iCell  ); j_index.push_back(jCell-2);
    i_index.push_back(iCell+1); j_index.push_back(jCell-2);
    i_index.push_back(iCell+2); j_index.push_back(jCell-2);
    i_index.push_back(iCell-2); j_index.push_back(jCell-1);
    i_index.push_back(iCell+2); j_index.push_back(jCell-1);
    i_index.push_back(iCell-2); j_index.push_back(jCell  );
    i_index.push_back(iCell+2); j_index.push_back(jCell  );
    i_index.push_back(iCell-2); j_index.push_back(jCell+1);
    i_index.push_back(iCell+2); j_index.push_back(jCell+1);
    i_index.push_back(iCell-2); j_index.push_back(jCell+2);
    i_index.push_back(iCell-1); j_index.push_back(jCell+2);
    i_index.push_back(iCell  ); j_index.push_back(jCell+2);
    i_index.push_back(iCell+1); j_index.push_back(jCell+2);
    i_index.push_back(iCell+2); j_index.push_back(jCell+2);

    /* Third ring incomplete (i.e. skip the corners) */
    i_index.push_back(iCell-2); j_index.push_back(jCell-3);
    i_index.push_back(iCell-1); j_index.push_back(jCell-3);
    i_index.push_back(iCell  ); j_index.push_back(jCell-3);
    i_index.push_back(iCell+1); j_index.push_back(jCell-3);
    i_index.push_back(iCell+2); j_index.push_back(jCell-3);
    i_index.push_back(iCell-3); j_index.push_back(jCell-2);
    i_index.push_back(iCell+3); j_index.push_back(jCell-2);
    i_index.push_back(iCell-3); j_index.push_back(jCell-1);
    i_index.push_back(iCell+3); j_index.push_back(jCell-1);
    i_index.push_back(iCell-3); j_index.push_back(jCell  );
    i_index.push_back(iCell+3); j_index.push_back(jCell  );
    i_index.push_back(iCell-3); j_index.push_back(jCell+1);
    i_index.push_back(iCell+3); j_index.push_back(jCell+1);
    i_index.push_back(iCell-3); j_index.push_back(jCell+2);
    i_index.push_back(iCell+3); j_index.push_back(jCell+2);
    i_index.push_back(iCell-2); j_index.push_back(jCell+3);
    i_index.push_back(iCell-1); j_index.push_back(jCell+3);
    i_index.push_back(iCell  ); j_index.push_back(jCell+3);
    i_index.push_back(iCell+1); j_index.push_back(jCell+3);
    i_index.push_back(iCell+2); j_index.push_back(jCell+3);

    break;

  case 1: // one ring of cells around (iCell,jCell)

    i_index.push_back(iCell);   j_index.push_back(jCell); /* cell (iCell,jCell) */

    /* First ring */
    i_index.push_back(iCell-1); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell-1);
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell+1); j_index.push_back(jCell+1);

    /* Second ring incomplete (i.e. skip the corners) */
    i_index.push_back(iCell-1); j_index.push_back(jCell-2);
    i_index.push_back(iCell  ); j_index.push_back(jCell-2);
    i_index.push_back(iCell+1); j_index.push_back(jCell-2);
    i_index.push_back(iCell-2); j_index.push_back(jCell-1);
    i_index.push_back(iCell+2); j_index.push_back(jCell-1);
    i_index.push_back(iCell-2); j_index.push_back(jCell  );
    i_index.push_back(iCell+2); j_index.push_back(jCell  );
    i_index.push_back(iCell-2); j_index.push_back(jCell+1);
    i_index.push_back(iCell+2); j_index.push_back(jCell+1);
    i_index.push_back(iCell-1); j_index.push_back(jCell+2);
    i_index.push_back(iCell  ); j_index.push_back(jCell+2);
    i_index.push_back(iCell+1); j_index.push_back(jCell+2);

    break;

  default: // general expression
    throw runtime_error("HighOrder2D<SOLN_STATE>::getEnlargedReconstructionStencil() doesn't support the current number of rings!");
  }//endswitch  
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the deviated (i.e. different than central) reconstruction stencil of 
 * cell (iCell,jCell).
 * Use the number of rings and the class variables caring information
 * about the opaqueness/transparency of boundaries to determine how far the stencil extends.
 * This routine DOES'T generate a central stencil if the cell gets affected by the presence of special boundaries!!!
 * It is also MORE EXPENSIVE than the regular routine that generates central stencils!!!
 * The stencil is biased to the mesh interior and it might be extended further than a central stencil.
 *
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 * \param [in] IsStencilExtended flag for extending or not the stencil. By default is true.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetDeviatedReconstructionStencil(const int &iCell, const int &jCell,
							       IndexType & i_index, IndexType & j_index,
							       const int &rings,
							       bool IsStencilExtended) const{

  // Reset indexes
  i_index.clear();
  j_index.clear();
  
  // The space around the cell (iCell,jCell) is divided into 8 regions relative to the cell faces (i.e. NW,W,SW,S,SE,E,NE,N).
  // There are 12 variable indexes which can fully control how the stencil is formed.

  int i,j;
  int Imin_NW, Jmax_NW;	// North-West region
  int Imin_W;	        // West region
  int Imin_SW, Jmin_SW;	// South-West region
  int Jmin_S;	        // South region
  int Imax_SE, Jmin_SE;	// South-East region
  int Imax_E;	        // East region
  int Imax_NE, Jmax_NE;	// North-East region
  int Jmax_N;	        // North region

  /* Set indexes as if a central stencil can be set. */
  Imin_NW = iCell-rings;
  Jmax_NW = jCell+rings;
  Imin_W  = iCell-rings;
  Imin_SW = iCell-rings;
  Jmin_SW = jCell-rings;
  Jmin_S  = jCell-rings;
  Imax_SE = iCell+rings;
  Jmin_SE = jCell-rings;
  Imax_E  = iCell+rings;
  Imax_NE = iCell+rings;
  Jmax_NE = jCell+rings;
  Jmax_N  = jCell+rings;

  
  // Decide whether the stencil is extended or not.
  // An extended stencil will add cells in the opposite direction of the restricted boundary.
  // By default the stencil is extended.
  if ( (CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS == ON &&
	CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL == OFF) ){
    
    IsStencilExtended = false;
  }


  /* ===  Cell to the left of West block boundary === 
   * Obs: This cell is unrestricted to West.
   */
  if ( iCell < ICl ){

    // ==== Cover cells with iCell < ICl ====
      
    if ( jCell < JCl){                // Case A
      // South Bnd influence
      if (W_SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl - 1){
	  Jmax_NW = Jmax_N = JCl-1;     /* limit Jmax */

	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    /* ensure valid Jmin */
	    Jmin_SW = max(Jmin_SW, 0);
	    Jmin_S  = max(Jmin_S, 0);
	  }
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCl-1;     /* limit Jmax */

	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;
	    /* ensure valid Jmin */
	    Jmin_SW = max(Jmin_SW, 0);
	    Jmin_S  = max(Jmin_S, 0); 
	    Jmin_SE = max(Jmin_SE, 0);
	  }
	}
      } else if (S_WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl - 1){
	  Imax_E = Imax_SE = ICl-1;     /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	    /* ensure valid Imin */
	    Imin_W  = max(Imin_W, 0);
	    Imin_SW = max(Imin_SW, 0);
	  }
	} else {
	  Imax_E = Imax_SE = Imax_NE = ICl-1;     /* limit Imax */
	    
	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	    Imin_NW -= 1;
	    /* ensure valid Imin */
	    Imin_W  = max(Imin_W, 0);
	    Imin_SW = max(Imin_SW, 0);
	    Imin_NW = max(Imin_NW, 0);
	  }
	}
      }

    } else if (jCell < JCl+rings){    // Case B
      // South Bnd influence
      if (W_SouthBnd.IsReconstructionStencilAffected() && SouthBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;         /* limit Jmin */
	  
	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1;
	}
      } else if (W_SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl - 1){
	  Jmin_SW = Jmin_S = JCl;                 /* limit Jmin */
	    
	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	  }
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCl;         /* limit Jmin */
	    
	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	  }
	}
      } else if (SouthBnd.IsReconstructionStencilAffected()){
	if (ICl-1-iCell > jCell - JCl){
	  Imax_SE = ICl - 1;                       /* limit Imax */

	  /* DON'T extend Imin due to smoothness indicator calculation reasons */
	} else {
	  Jmin_SE = JCl;                            /* limit Jmin */

	  /* DON'T extend Jmax due to smoothness indicator calculation reasons */
	}
      }
	
    } else if (jCell>JCu-rings && jCell <= JCu){	        // Case C (JCu-rings < jCell <= JCu)
      // North Bnd influence
      if (W_NorthBnd.IsReconstructionStencilAffected() && NorthBnd.IsReconstructionStencilAffected()){
	Jmax_NW = Jmax_N = Jmax_NE = JCu;         /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;
	}
      } else if (W_NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl -1){
	  Jmax_NW = Jmax_N = JCu;                 /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	  }
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCu;         /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1; 
	  } 
	}
      } else if (NorthBnd.IsReconstructionStencilAffected()){
	if (ICl-1-iCell > JCu-jCell){
	  Imax_NE = ICl - 1;                       /* limit Imax */
	    
	  /* DON'T extend Imin due to smoothness indicator calculation reasons */
	} else {
	  Jmax_NE = JCu;                            /* limit Jmax */

	  /* DON'T extend Jmin due to smoothness indicator calculation reasons */
	}
      }

    } else if (jCell > JCu){	                 // Case D (jCell > JCu)
      // North Bnd influence
      if (W_NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl - 1){
	  Jmin_SW = Jmin_S = JCu+1;     /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    /* ensure valid Jmax */
	    Jmax_NW = min(Jmax_NW, JCu+Ng);
	    Jmax_N  = min(Jmax_N, JCu+Ng);	    
	  }
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCu+1;     /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	    /* ensure valid Jmax */
	    Jmax_NW = min(Jmax_NW, JCu+Ng);
	    Jmax_N  = min(Jmax_N, JCu+Ng);
	    Jmax_NE = min(Jmax_NE, JCu+Ng);
	  }
	}
      } else if (N_WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu + 1){
	  Imax_NE = Imax_E = ICl-1;     /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    /* ensure valid Imin */
	    Imin_NW = max(Imin_NW, 0);
	    Imin_W  = max(Imin_W, 0);
	  }
	} else {
	  Imax_NE = Imax_E = Imax_SE = ICl-1;     /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	    /* ensure valid Imin */
	    Imin_NW = max(Imin_NW, 0);
	    Imin_W  = max(Imin_W, 0);
	    Imin_SW = max(Imin_SW, 0);
	  }
	}
      }
    }	// endif (Case D)

  } else if (iCell < ICl + rings) {

    // ==== Cover cells with (ICl <= iCell < ICl+rings) ====

    if ( jCell < JCl){	         // Case A
      // West Bnd influence
      if (S_WestBnd.IsReconstructionStencilAffected() && WestBnd.IsReconstructionStencilAffected()){
	Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */

	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_SE += 1;
	  Imax_E  += 1;
	  Imax_NE += 1;
	}
      } else if (S_WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl -1){
	  Imin_SW = Imin_W = ICl;     /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_SE += 1;
	    Imax_E  += 1;
	  }
	} else {
	  Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_SE += 1;
	    Imax_E  += 1;
	    Imax_NE += 1;  
	  }  
	}
      } else if (WestBnd.IsReconstructionStencilAffected()){
	if (iCell-ICl > JCl-1-jCell){
	  Imin_NW = ICl;              /* limit Imin */
	  
	  /* DON'T extend Imax due to smoothness indicator calculation reasons */
	} else {
	  // South Bnd influence
	  Jmax_NW = JCl - 1;                     /* limit Jmax */
	    
	  /* DON'T extend Jmin due to smoothness indicator calculation reasons */
	}
      }

    } else if (jCell < JCl + rings){   // Case B
      if (WestBnd.IsReconstructionStencilAffected() && SouthBnd.IsReconstructionStencilAffected()){
	Imin_SW = Imin_W = Imin_NW = ICl;      /* limit Imin */
	Jmin_SW = Jmin_S = Jmin_SE = JCl;      /* limit Jmin */

	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_SE += 1;
	  Imax_E  += 1;
	  Imax_NE += 1;
	  
	  /* extend Jmax */
	  Jmax_NW += 1; 
	  Jmax_N  += 1;
	  Jmax_NE += 1;
	}
      } else if (WestBnd.IsReconstructionStencilAffected() && S_WestBnd.IsReconstructionStencilAffected()){
	Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */
  
	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_SE += 1;
	  Imax_E  += 1;
	  Imax_NE += 1;	
	}  
      } else if (SouthBnd.IsReconstructionStencilAffected() && W_SouthBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;     /* limit Jmin */

	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1;
	}
      } else if (WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl){
	  Imin_W = Imin_NW = ICl;             /* limit Imin */

	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_NE += 1;
	    Imax_E  += 1;
	  }
	} else {
	  Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_SE += 1;
	    Imax_E  += 1;
	    Imax_NE += 1;
	  } 
	}
      } else if (SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl){
	  Jmin_S = Jmin_SE = JCl;              /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	  }
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCl;              /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	  }
	}
      } else if (W_SouthBnd.IsReconstructionStencilAffected()){ // S_WestBnd is also affecting the stencil!
	if (iCell > jCell){
	  // Limit to West side
	  Imin_SW = ICl;
	  // DON'T extend stencil
	} else {
	  // Limit to South side
	  Jmin_SW = JCl;
	  // DON'T extend stencil
	}
      }

    } else if (jCell <= JCu - rings){  // Case C 
      // West Bnd influence
      if (WestBnd.IsReconstructionStencilAffected()){
	Imin_NW = Imin_W = Imin_SW = ICl;              /* limit Imin */

	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_NE += 1;
	  Imax_E  += 1;
	  Imax_SE += 1;
	}
      }

    } else if (jCell <= JCu){	         // Case D
      if (WestBnd.IsReconstructionStencilAffected() && NorthBnd.IsReconstructionStencilAffected()){
	Imin_NW = Imin_W = Imin_SW = ICl;             /* limit Imin */
	Jmax_NW = Jmax_N = Jmax_NE = JCu;             /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_NE += 1;
	  Imax_E  += 1;
	  Imax_SE += 1;
	  
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;
	}

      } else if (WestBnd.IsReconstructionStencilAffected() && N_WestBnd.IsReconstructionStencilAffected()){
	Imin_NW = Imin_W = Imin_SW = ICl;             /* limit Imin */

	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_NE += 1;
	  Imax_E  += 1;
	  Imax_SE += 1;	
	}  
      } else if (NorthBnd.IsReconstructionStencilAffected() && W_NorthBnd.IsReconstructionStencilAffected()){
	Jmax_NW = Jmax_N = Jmax_NE = JCu;             /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;
	}
      } else if (WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu){
	  Imin_W = Imin_SW = ICl;             /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_E  += 1;
	    Imax_SE += 1;	 
	  } 	    
	} else {
	  Imin_NW = Imin_W = Imin_SW = ICl;             /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_NE += 1;
	    Imax_E  += 1;
	    Imax_SE += 1;	 
	  } 	    
	}
      } else if (NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICl){
	  Jmax_N = Jmax_NE = JCu;             /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;	 
	  }   	    
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCu;             /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;	 
	  }   
	}
      } else if (W_NorthBnd.IsReconstructionStencilAffected()){ // N_WestBnd also affects the stencil
	if (iCell-ICl > JCu-jCell){
	  // Limit to West side
	  Imin_NW = ICl;
	  // DON'T extend stencil
	} else {
	  // Limit to North side
	  Jmax_NW = JCu;
	  // DON'T extend stencil
	}
      }

    } else {			         // Case E (jCell > JCu)
      // West Bnd influence
      if (N_WestBnd.IsReconstructionStencilAffected() && WestBnd.IsReconstructionStencilAffected()){
	Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */
	  
	if (IsStencilExtended){
	  /* extend Imax */
	  Imax_SE += 1;
	  Imax_E  += 1;
	  Imax_NE += 1;	
	}  
      } else if (N_WestBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu + 1){
	  Imin_W = Imin_NW = ICl;     /* limit Imin */
	    
	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_E  += 1;
	    Imax_NE += 1;	 
	  } 	    
	} else {
	  Imin_SW = Imin_W = Imin_NW = ICl;     /* limit Imin */

	  if (IsStencilExtended){	    
	    /* extend Imax */
	    Imax_SE += 1;
	    Imax_E  += 1;
	    Imax_NE += 1;	 
	  } 
	}
      } else if (WestBnd.IsReconstructionStencilAffected()){
	if (iCell-ICl > jCell-(JCu+1)){	    
	  Imin_SW = ICl;                         /* limit Imin */

	  /* DON'T extend Imax due to smoothness indicator calculation reasons */
	} else {
	  Jmin_SW = JCu + 1;                     /* limit Jmin */

	  /* DON'T extend Jmax due to smoothness indicator calculation reasons */
	}
      }
    }	// endif (Case E)

  } else if (iCell <= ICu - rings) {

    // ==== Cover cells with (ICl+rings <= iCell <= ICu-rings) ====

    if ( jCell < JCl + rings){         // Case A
      // South Bnd influence
      if (SouthBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;     /* limit Jmin */

	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1;
	}
      }

    } else if (jCell > JCu - rings){   // Case B
      // North Bnd influence
      if (NorthBnd.IsReconstructionStencilAffected()){
	Jmax_NW = Jmax_N = Jmax_NE = JCu;     /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;	
	}  
      }
    }	// endif (Case B)

  } else if (iCell <= ICu) {

    // ==== Cover cells with (ICu-rings < iCell <= ICu) ====

    if ( jCell < JCl){	         // Case A
      // East Bnd influence
      if (S_EastBnd.IsReconstructionStencilAffected() && EastBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;     /* limit Imax */

	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;
	}
      } else if (S_EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl-1){
	  Imax_E = Imax_SE = ICu;     /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_W  -= 1;
	    Imin_SW -= 1;  
	  } 
	} else {
	  Imax_NE = Imax_E = Imax_SE = ICu;     /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	  }
	}
      } else if (EastBnd.IsReconstructionStencilAffected()){
	if (ICu-iCell > JCl-1-jCell){
	  Imax_NE = ICu;                         /* limit Imax */

	  /* DON'T extend Imin due to smoothness indicator calculation reasons */
	} else {
	  // South Bnd influence
	  Jmax_NE = JCl - 1;                     /* limit Jmax */

	  /* DON'T extend Jmin due to smoothness indicator calculation reasons */	  
	}
      }

    } else if (jCell < JCl + rings){   // Case B
      if (SouthBnd.IsReconstructionStencilAffected() && EastBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;        /* limit Jmin */
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1;

	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;
	}
      } else if (EastBnd.IsReconstructionStencilAffected()  && S_EastBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;
	}
      } else if (SouthBnd.IsReconstructionStencilAffected() && E_SouthBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;        /* limit Jmin */

	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1; 
	}
      } else if (EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl){
	  Imax_NE = Imax_E = ICu;              /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	  }
	} else {
	  Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	  }
	}
      } else if (SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu){
	  Jmin_SW = Jmin_S = JCl;     /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	  }
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCl;    /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;  
	  }   
	}
      } else if (E_SouthBnd.IsReconstructionStencilAffected()){ // S_EastBnd also affects the stencil
	if ( (ICu-iCell) > (jCell - JCl) ){
	  // Limit to East side
	  Imax_SE = ICu;
	  // DON'T extend stencil
	} else {
	  // Limit to South side 
	  Jmin_SE = JCl;
	  // DON'T extend stencil
	}
      }

    } else if (jCell <= JCu - rings){  // Case C 
      if (EastBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;
	}
      }

    } else if (jCell <= JCu){	         // Case D
      if (EastBnd.IsReconstructionStencilAffected() && NorthBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */
	Jmax_NW = Jmax_N = Jmax_NE = JCu;        /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;	  

	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;
	}

      } else if (EastBnd.IsReconstructionStencilAffected()  && N_EastBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;
	}
      } else if (NorthBnd.IsReconstructionStencilAffected() && E_NorthBnd.IsReconstructionStencilAffected()){
	Jmax_NW = Jmax_N = Jmax_NE = JCu;        /* limit Jmax */

	if (IsStencilExtended){
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;
	}
      } else if (EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu){
	  Imax_E = Imax_SE = ICu;             /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	  }
	} else {
	  Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */

	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	  }
	}
      } else if (NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu){
	  Jmax_NW = Jmax_N = JCu;        /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	  }
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCu;        /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;
	  }
	}
      } else if (E_NorthBnd.IsReconstructionStencilAffected()){ // N_EastBnd also affects the stencil
	if (ICu-iCell > JCu-jCell){
	  // Limit to East side
	  Imax_NE = ICu;
	  // DON'T extend stencil
	} else {
	  // Limit to North side
	  Jmax_NE = JCu;
	  // DON'T extend stencil
	}
      }

    } else {			         // Case E (jCell > JCu)
      // East Bnd influence
      if (N_EastBnd.IsReconstructionStencilAffected() && EastBnd.IsReconstructionStencilAffected()){
	Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */
	  
	if (IsStencilExtended){
	  /* extend Imin */
	  Imin_NW -= 1;
	  Imin_W  -= 1;
	  Imin_SW -= 1;	
	}  
      } else if (N_EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu + 1){
	  Imax_NE = Imax_E = ICu;        /* limit Imax */
	    
	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;  
	  }  
	} else {
	  Imax_NE = Imax_E = Imax_SE = ICu;        /* limit Imax */
	    
	  if (IsStencilExtended){
	    /* extend Imin */
	    Imin_NW -= 1;
	    Imin_W  -= 1;
	    Imin_SW -= 1;
	  }
	}
      } else if (EastBnd.IsReconstructionStencilAffected()){
	if (ICu - iCell > jCell - (JCu+1)) {
	  Imax_SE = ICu;                         /* limit Imax */
	    
	  /* DON'T extend Imin due to smoothness indicator calculation reasons */
	} else {
	  // North Bnd influence
	  Jmin_SE = JCu + 1;                     /* limit Jmin */
	    
	  /* DON'T extend Jmax due to smoothness indicator calculation reasons */
	}
      }
    }// endif (Case E)

  } else {

    // ==== Cover cells with (iCell > ICu) ====
    if ( jCell < JCl){	         // Case A
      // South Bnd influence
      if (E_SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu + 1){
	  Jmax_N = Jmax_NE = JCl-1;     /* limit Jmax */

	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;
	    /* ensure valid Jmin */
	    Jmin_S  = max(Jmin_S, 0);
	    Jmin_SE = max(Jmin_SE, 0);
	  }
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCl-1;     /* limit Jmax */

	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;
	    /* ensure valid Jmin */
	    Jmin_SW = max(Jmin_SW, 0);
	    Jmin_S  = max(Jmin_S, 0);
	    Jmin_SE = max(Jmin_SE, 0);
	  }
	}
      } else if (S_EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCl - 1){
	  Imin_W = Imin_SW = ICu+1;              /* limit Imin */

	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_E  += 1;
	    Imax_SE += 1;
	    /* ensure valid Imax */
	    Imax_E  = min(Imax_E, ICu+Ng);
	    Imax_SE = min(Imax_SE, ICu+Ng);	    
	  }
	} else {
	  Imin_NW = Imin_W = Imin_SW = ICu+1;    /* limit Imin */

	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_NE += 1;
	    Imax_E  += 1;
	    Imax_SE += 1;
	    /* ensure valid Imax */
	    Imax_NE = min(Imax_NE, ICu+Ng);	    
	    Imax_E  = min(Imax_E, ICu+Ng);
	    Imax_SE = min(Imax_SE, ICu+Ng);	    
	  }
	}
      }

    } else if (jCell < JCl+rings){     // Case B
      // South Bnd influence
      if (E_SouthBnd.IsReconstructionStencilAffected() && SouthBnd.IsReconstructionStencilAffected()){
	Jmin_SW = Jmin_S = Jmin_SE = JCl;         /* limit Jmin */
	  
	if (IsStencilExtended){
	  /* extend Jmax */
	  Jmax_NW += 1;
	  Jmax_N  += 1;
	  Jmax_NE += 1;	
	}  
      } else if (E_SouthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu + 1){
	  Jmin_S = Jmin_SE = JCl;         /* limit Jmin */
	    
	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_N  += 1;
	    Jmax_NE += 1;	 
	  } 
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCl;         /* limit Jmin */
	    
	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;	 
	  } 
	}
      } else if (SouthBnd.IsReconstructionStencilAffected()){
	if (iCell - (ICu+1) > jCell-JCl) {
	  // East Bnd influence
	  Imin_SW = ICu + 1;                       /* limit Imin */
	    
	  /* DON'T extend Imax due to smoothness indicator calculation reasons */
	} else {
	  Jmin_SW = JCl;                            /* limit Jmin */
	    
	  /* DON'T extend Jmax due to smoothness indicator calculation reasons */	  
	}
      } // endif (Case B)

    } else if (jCell>JCu-rings && jCell <= JCu){	        // Case C (JCu-rings < jCell <= JCu)
      // North Bnd influence
      if (E_NorthBnd.IsReconstructionStencilAffected() && NorthBnd.IsReconstructionStencilAffected()){
	Jmax_NW = Jmax_N = Jmax_NE = JCu;         /* limit Jmax */
	  
	if (IsStencilExtended){
	  /* extend Jmin */
	  Jmin_SW -= 1;
	  Jmin_S  -= 1;
	  Jmin_SE -= 1;	
	}  
      } else if (E_NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu + 1){
	  Jmax_N = Jmax_NE = JCu;         /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;	 
	  } 
	} else {
	  Jmax_NW = Jmax_N = Jmax_NE = JCu;         /* limit Jmax */
	    
	  if (IsStencilExtended){
	    /* extend Jmin */
	    Jmin_SW -= 1;
	    Jmin_S  -= 1;
	    Jmin_SE -= 1;	 
	  } 
	}
      } else if (NorthBnd.IsReconstructionStencilAffected()){
	if (iCell - (ICu+1) > JCu - jCell) {
	  // East Bnd influence
	  Imin_NW = ICu + 1;                        /* limit Imax */
	    
	  /* DON'T extend Imin due to smoothness indicator calculation reasons */
	} else {
	  Jmax_NW = JCu;                         /* limit Jmax */
	    
	  /* DON'T extend Jmin due to smoothness indicator calculation reasons */ 
	}
      } // endif (Case C)

    } else if (jCell > JCu){	                        // Case D (jCell > JCu)
      // North Bnd influence
      if (E_NorthBnd.IsReconstructionStencilAffected()){
	if (iCell == ICu + 1){
	  Jmin_S = Jmin_SE = JCu+1;     /* limit Jmin */

	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	    /* ensure valid Jmax */
	    Jmax_N  = min(Jmax_N, JCu+Ng);
	    Jmax_NE = min(Jmax_NE, JCu+Ng);	    
	  }
	} else {
	  Jmin_SW = Jmin_S = Jmin_SE = JCu+1;     /* limit Jmin */
	    
	  if (IsStencilExtended){
	    /* extend Jmax */
	    Jmax_NW += 1;
	    Jmax_N  += 1;
	    Jmax_NE += 1;
	    /* ensure valid Jmax */
	    Jmax_NW = min(Jmax_NW, JCu+Ng);	    
	    Jmax_N  = min(Jmax_N, JCu+Ng);
	    Jmax_NE = min(Jmax_NE, JCu+Ng);	    
	  }
	}
      } else if (N_EastBnd.IsReconstructionStencilAffected()){
	if (jCell == JCu + 1){
	  Imin_NW = Imin_W = ICu + 1;   /* limit Imin */

	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_NE += 1;
	    Imax_E  += 1;
	    /* ensure valid Imax */
	    Imax_NE = min(Imax_NE, ICu+Ng);
	    Imax_E  = min(Imax_E, ICu+Ng);
	  }
	} else {
	  Imin_NW = Imin_W = Imin_SW = ICu + 1;   /* limit Imin */

	  if (IsStencilExtended){
	    /* extend Imax */
	    Imax_NE += 1;
	    Imax_E  += 1;
	    Imax_SE += 1;
	    /* ensure valid Imax */
	    Imax_NE = min(Imax_NE, ICu+Ng);
	    Imax_E  = min(Imax_E, ICu+Ng);
	    Imax_SE = min(Imax_SE, ICu+Ng);
	  }
	}
      }
    }	// endif (Case D)

  } // endif (iCell)

  /* ===== Form stencil ===== */
  i_index.push_back(iCell);
  j_index.push_back(jCell);

  // Add cells from NW region
  for (i=Imin_NW; i<=iCell-1; ++i){
    for (j=jCell+1; j<=Jmax_NW; ++j){
      i_index.push_back(i);
      j_index.push_back(j);
    }// endif
  }// endif

  // Add cells from W region
  for (i=Imin_W; i<=iCell-1; ++i){
    i_index.push_back(i);
    j_index.push_back(jCell);
  }// endif

  // Add cells from SW region
  for (i=Imin_SW; i<=iCell-1; ++i){
    for (j=Jmin_SW; j<=jCell-1; ++j){
      i_index.push_back(i);
      j_index.push_back(j);
    }// endif
  }// endif

  // Add cells from S region
  for (j=Jmin_S; j<=jCell-1; ++j){
    i_index.push_back(iCell);
    j_index.push_back(j);
  }// endif

  // Add cells from SE region
  for (i=iCell+1; i<=Imax_SE; ++i){
    for (j=Jmin_SE; j<=jCell-1; ++j){
      i_index.push_back(i);
      j_index.push_back(j);
    }// endif
  }// endif

  // Add cells from E region
  for (i=iCell+1; i<=Imax_E; ++i){
    i_index.push_back(i);
    j_index.push_back(jCell);
  }// endif

  // Add cells from NE region
  for (i=iCell+1; i<=Imax_NE; ++i){
    for (j=jCell+1; j<=Jmax_NE; ++j){
      i_index.push_back(i);
      j_index.push_back(j);
    }// endif
  }// endif

  // Add cells from N region
  for (j=jCell+1; j<=Jmax_N; ++j){
    i_index.push_back(iCell);
    j_index.push_back(j);
  }// endif

}

template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::displayDeviatedReconstructionStencil(ostream & out,
								   const int &iCell, const int &jCell,
								   const int &rings) const{

  // Set local variables
  int i,j;
  IndexType i_index, j_index;
  char** StencilMap;

  // Allocate StencilMap memory
  StencilMap = new char* [9];
  for (i = 0; i<9; ++i){
    StencilMap[i] = new char [9];
    // Initialize StencilMap to emply character
    for (j = 0; j<9; ++j){
      StencilMap[i][j] = ' ';
    }
  }

  // Mark the reconstructed cell with the '@' sign
  StencilMap[4][4] = '@';

  // Determine the deviated stencil for cell (iCell,jCell)
  SetDeviatedReconstructionStencil(iCell, jCell, i_index, j_index, rings);  

  // Mark the presence of each neighbour cell in the stencil with '*' sign
  for (i = 1; i<i_index.size(); ++i){
    StencilMap[4 + (i_index[i] - iCell)][4 + (j_index[i] - jCell)] = '*';
  }

  // Output StencilMap
  for (j = 8; j>=0; --j){
    for (i = 0; i<=8; ++i){
      out << StencilMap[i][j];
    }
    out << endl;
  }

}

/*!
 * Return the type of reconstruction for a given cell (iCell,jCell).
 */
template<class SOLN_STATE> inline
char HighOrder2D<SOLN_STATE>::getCellReconstructionType(const int& iCell, const int& jCell) const {
  if (ReconstructionTypeMap != NULL){
    return ReconstructionTypeMap[iCell][jCell];
  } else {
    return 'n';
  }
}

/*! 
 * Compute the gradient at an inter-cellular face using the first 
 * solution parameter based on the Green-Gauss reconstruction.
 * The path used to compute the gradient is formed by the centroids 
 * of the first order neighbours.
 * As an example, for EAST face the cells are: 0, S, SE, E, NE, N, 0
 * where 0 is the (iCell,jCell) cell and the orientation of the other 
 * cells is relative to 0 cell. 
 * The last index repeats the first one to show that the path is closed. 
 *
 * \param [in] Face Specify which face the gradient is computed at.
 * \param [out] GradU_face the gradient value is returned here 
 *
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::
GreenGauss_FaceGradient_CentroidPathCartesianMesh(Soln_Block_Type &SolnBlk,
						  const int &iCell, const int &jCell,
						  const int &Face,
						  Vector2D & GradU_face,
						  const Soln_State &
						  (Soln_Block_Type::*ReconstructedSoln)(const int &,const int &) const){

  int StencilSize(7);		// the last entry is to close the path

  int n;
  int Info;
  double PolygonArea, Length;
  Vector2D Centroids[StencilSize], PolygonCentroid;
  Vector2D normal;
  IndexType i_index, j_index;
  double Um;		// average solution along the current integration segment

  i_index.reserve(StencilSize);
  j_index.reserve(StencilSize);

  // Form the supporting stencil differently for each face
  switch (Face){

  case NORTH:
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell-1); j_index.push_back(jCell+1);
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    break;

  case SOUTH:
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    break;

  case WEST:
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell-1); j_index.push_back(jCell+1);
    i_index.push_back(iCell-1); j_index.push_back(jCell  );
    i_index.push_back(iCell-1); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell  ); j_index.push_back(jCell  );

    break;

  case EAST:
    i_index.push_back(iCell  ); j_index.push_back(jCell  );
    i_index.push_back(iCell  ); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell-1);
    i_index.push_back(iCell+1); j_index.push_back(jCell  );
    i_index.push_back(iCell+1); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell+1);
    i_index.push_back(iCell  ); j_index.push_back(jCell  );

    break;
  }


  // Form the array of the centroids
  for (n = 0; n < i_index.size(); ++n){
    Centroids[n] = CellCenter(i_index[n], j_index[n]);
  }

  // Determine the area and centroid of the closed centroid path
  Info = polyCentroid(Centroids, StencilSize, PolygonCentroid, PolygonArea);

  // Reset gradient
  GradU_face = Vector2D(0);

  // Calculate average gradient
  for (n = 1; n < i_index.size(); ++n){

    // Calculate segment length between centroids
    Length = abs(Centroids[n] - Centroids[n-1]);

    // Calculate normal
    normal = Vector2D( (Centroids[n].y - Centroids[n-1].y),
		       -(Centroids[n].x - Centroids[n-1].x ) )/Length;

    // Calculate average solution with the first parameter
    Um = 0.5*( (SolnBlk.*ReconstructedSoln)(i_index[n  ], j_index[n  ])[1] +
	       (SolnBlk.*ReconstructedSoln)(i_index[n-1], j_index[n-1])[1] );

    // Add segment contribution
    GradU_face += normal * Um * Length;

  }

  GradU_face /= PolygonArea;

}

/*! 
 * Compute the integral over the block geometry of the error between the
 * reconstructed polynomial and the function provided as input. 
 *
 * \param [in] FuncObj  The function relative to which the error is evaluated.
 *                      It is assumed that the exact solution can take two arguments
 *                      (x & y position) and returns a double.
 * \param [in] parameter The state variable which is used for computing the errors
 * \param [in] digits  The targeted number of exact digits with which the integral is evaluated
 */
template<class SOLN_STATE> 
template<typename Function_Object_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeSolutionErrors(const Function_Object_Type FuncObj,
						    const unsigned &parameter,
						    const int &digits){
  
  int StartI_Int, EndI_Int, StartJ_Int, EndJ_Int; // Integration indexes for cells with straight edges.
  int i,j;
  bool _integrate_with_curved_boundaries(false);
  double ErrorTemp;
  char Case;
  
  /* Algorithm:
     1. Decide where integration for straight quads can be used.
        If low boundary representation is used than all quads are straight.
	If high-order boundary representation is provided, then all cells near boundaries are 
        considered to be curved, even if they might actually have boundary straight lines.
     2. Integrate over straight quads.
     3. Integrate over the curved cells if necessary.
  */

  /* Possible cases:
     
     a. Error based on entropy:  parameter = 0;
     b. Error based on a particular parameter:  parameter = 1-NumberOfVariables();
   */

  if (parameter == 0){
    Case = 'a';
  } else if (parameter >= 1 && parameter <= NumberOfVariables() ){
    Case = 'b';
  }

  // Decide the range of integration
  if (Geom->IsHighOrderBoundary()){
    if (CENO_Execution_Mode::IGNORE_CURVED_BOUNDARIES_FOR_ACCURACY_ASSESSMENT){
      // Check if West spline has more than 2 control points
      if (Geom->BndWestSpline.np > 2){
	// Ignore cells near this boundary
	StartI_Int = ICl + 1;
      } else {
	// Include cells near this boundary
	StartI_Int = ICl;
      }

      // Check if East spline has more than 2 control points
      if (Geom->BndEastSpline.np > 2){
	// Ignore cells near this boundary
	EndI_Int   = ICu - 1;
      } else {
	// Include cells near this boundary
	EndI_Int   = ICu;
      }

      // Check if South spline has more than 2 control points
      if (Geom->BndSouthSpline.np > 2){
	// Ignore cells near this boundary
	StartJ_Int = JCl + 1;
      } else {
	// Include cells near this boundary
	StartJ_Int = JCl;
      }

      // Check if North spline has more than 2 control points
      if (Geom->BndNorthSpline.np > 2){
	// Ignore cells near this boundary
	EndJ_Int   = JCu - 1;
      } else {
	// Include cells near this boundary
	EndJ_Int   = JCu;
      }

    } else {
      _integrate_with_curved_boundaries = true;
      throw runtime_error("HighOrder2D<SOLN_STATE>::ComputeSolutionErrors() Warning! Integration with curved bnds is not setup!");
    }
  } else {
    StartI_Int = ICl;
    EndI_Int   = ICu;
    StartJ_Int = JCl;
    EndJ_Int   = JCu;
  }

  // Reset the error values
  ResetErrors();

  // Sum up the contribution from each straight quad cell
  for (j = StartJ_Int; j <= EndJ_Int; ++j){
    for (i = StartI_Int; i <= EndI_Int; ++i, ++CellsUsed){

      switch (Case){
      case 'a': 
	ErrorTemp = ComputeSolutionEntropyErrorL1(i,j,
						  FuncObj,
						  digits);

	ErrorL2 += ComputeSolutionEntropyErrorL2(i,j,
						 FuncObj,
						 digits);
	break;
	
      case 'b':
	ErrorTemp = ComputeSolutionErrorL1(i,j,
					   FuncObj,
					   parameter,
					   digits);

	ErrorL2 += ComputeSolutionErrorL2(i,j,
					  FuncObj,
					  parameter,
					  digits);
	break;
      }

      ErrorL1 += ErrorTemp;
      ErrorMax = max(ErrorMax, ErrorTemp/Geom->CellArea(i,j));

      TotalBlockArea += Geom->CellArea(i,j);
    }// endfor
  }// endfor

  // Sum up the contribution from the cells with curved boundaries
  if (_integrate_with_curved_boundaries){

    for (j = JCl; j <= JCu; ++j){
      // not available yet

      // West boundary
      // (ICl,j)

      // East boundary
      // (ICu,j)

      CellsUsed += 2;
    }// endfor

    for (i = ICl+1; i <= ICu-1; ++i){
      // not available yet

      // South boundary
      // (i,JCl)

      // North boundary
      // (i,JCu)

      CellsUsed += 2;
    }// endfor

  } //endif

}

/*! 
 * Compute the integral over the block geometry of the polynomial
 * reconstructions of two high-order variables over the whole domain. 
 *
 * \param [in] FuncObj  The function relative to which the error is evaluated
 * \param [in] parameter The parameter for which the reconstruction is evaluated
 * \param [in] digits  The targeted number of exact digits with which the integral is evaluated
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ComputeReconstructionsErrors(const HighOrder2D<Soln_State> & Obj,
							   const unsigned &parameter,
							   const int &digits){
  
  int StartI_Int, EndI_Int, StartJ_Int, EndJ_Int; // Integration indexes for cells with straight edges.
  int i,j;
  bool _integrate_with_curved_boundaries(false);
  double ErrorTemp;
  
  /* Algorithm:
     1. Decide where integration for straight quads can be used.
        If low boundary representation is used than all quads are straight.
	If high-order boundary representation is provided, then all cells near boundaries are 
        considered to be curved, even if they might actually have boundary straight lines.
     2. Integrate over straight quads.
     3. Integrate over the curved cells if necessary.
  */

  // Decide the range of integration
  if (Geom->IsHighOrderBoundary()){
    StartI_Int = ICl + 1;
    EndI_Int   = ICu - 1;
    StartJ_Int = JCl + 1;
    EndJ_Int   = JCu - 1;
    _integrate_with_curved_boundaries = true;
    throw runtime_error("HighOrder2D<SOLN_STATE>::ComputeReconstructionsErrors() Warning! Integration with curved bnds is not setup!");
  } else {
    StartI_Int = ICl;
    EndI_Int   = ICu;
    StartJ_Int = JCl;
    EndJ_Int   = JCu;
  }

  // Reset the error values
  ResetErrors();

  // Sum up the contribution from each straight quad cell
  for (j = StartJ_Int; j <= EndJ_Int; ++j){
    for (i = StartI_Int; i <= EndI_Int; ++i, ++CellsUsed){
      ErrorTemp = ComputeReconstructionsErrorL1(i,j,
					       Obj,
					       parameter,
					       digits);
 
      ErrorL1 += ErrorTemp;

      ErrorL2 += ComputeReconstructionsErrorL2(i,j,
					       Obj,
					       parameter,
					       digits);

      ErrorMax = max(ErrorMax, ErrorTemp/Geom->CellArea(i,j));

      TotalBlockArea += Geom->CellArea(i,j);
    }// endfor
  }// endfor

  // Sum up the contribution from the cells with curved boundaries
  if (_integrate_with_curved_boundaries){

    for (j = JCl; j <= JCu; ++j){
      // not available yet

      // West boundary
      // (ICl,j)

      // East boundary
      // (ICu,j)

      CellsUsed += 2;
    }// endfor

    for (i = ICl+1; i <= ICu-1; ++i){
      // not available yet

      // South boundary
      // (i,JCl)

      // North boundary
      // (i,JCu)

      CellsUsed += 2;
    }// endfor

  } //endif

}

/*! 
 * Compute the integral of the error function between what 
 * is considered as exact solution and the reconstructed 
 * polynomial over the domain of cell (iCell,jCell).
 * This result is used in the evaluation of L1 error norm.
 */
template<class SOLN_STATE>
template<typename Function_Object_Type> inline
double HighOrder2D<SOLN_STATE>::ComputeSolutionErrorL1(const int &iCell, const int &jCell,
						       const Function_Object_Type FuncObj,
						       const unsigned &parameter, const int & digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration
  Vector2D _dummy_Position(0.0); // defined only to provide the type of the position vector where the function is evaluated.

  return IntegrateOverTheCell(iCell,jCell,
			      error_function(FuncObj,
					     wrapped_member_function_one_parameter(this,
										   &ClassType::SolutionAtLocation,
										   _dummy_Position,
										   iCell, jCell,
										   parameter,
										   _dummy_result),
					     _dummy_result),
			      digits, _dummy_result);

}

/*! 
 * Compute the integral of the entropy error function between
 * a reference solution and the reconstructed 
 * polynomial over the domain of cell (iCell,jCell).
 * This result is used in the evaluation of L1 error norm.
 */
template<class SOLN_STATE>
template<typename Function_Object_Type> inline
double HighOrder2D<SOLN_STATE>::ComputeSolutionEntropyErrorL1(const int &iCell, const int &jCell,
							      const Function_Object_Type FuncObj,
							      const int & digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration

  return IntegrateOverTheCell(iCell,jCell,
			      error_function(FuncObj,
					     wrapped_soln_block_member_function(this,
										&ClassType::SolutionEntropyAtCoordinates,
										iCell, jCell,
										_dummy_result),
					     _dummy_result),
			      digits, _dummy_result);

}

/*! 
 * Compute the integral of the error function between the
 * polynomial reconstructions of two high-order variables
 * over the domain of cell (iCell,jCell).
 *
 * \note This subroutine doesn't check that the domains are the same! It's the caller's responsibility.
 */
template<class SOLN_STATE> inline
double HighOrder2D<SOLN_STATE>::ComputeReconstructionsErrorL1(const int &iCell, const int &jCell,
							      const HighOrder2D<Soln_State> & Obj,
							      const unsigned &parameter,
							      const int &digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration
  Vector2D _dummy_Position(0.0); // defined only to provide the type of the position vector where the function is evaluated.

  return IntegrateOverTheCell(iCell,jCell,
			      error_function(wrapped_member_function_one_parameter(&Obj,
										   &ClassType::SolutionAtLocation,
										   _dummy_Position,
										   iCell, jCell,
										   parameter,
										   _dummy_result),
					     wrapped_member_function_one_parameter(this,
										   &ClassType::SolutionAtLocation,
										   _dummy_Position,
										   iCell, jCell,
										   parameter,
										   _dummy_result),
					     _dummy_result),
			      digits, _dummy_result);
}

/*! 
 * Compute the integral of the squared error function between what 
 * is considered as exact solution and the reconstructed 
 * polynomial over the domain of cell (iCell,jCell).
 * This result is used in the evaluation of L2 error norm.
 */
template<class SOLN_STATE> 
template<typename Function_Object_Type> inline
double HighOrder2D<SOLN_STATE>::ComputeSolutionErrorL2(const int &iCell, const int &jCell,
						       Function_Object_Type FuncObj,
						       const unsigned &parameter, const int &digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration
  Vector2D _dummy_Position(0.0); // defined only to provide the type of the position vector where the function is evaluated.
  
  return IntegrateOverTheCell(iCell,jCell,
			      square_error_function(FuncObj,
						    wrapped_member_function_one_parameter(this,
											  &ClassType::SolutionAtLocation,
											  _dummy_Position,
											  iCell, jCell,
											  parameter,
											  _dummy_result),
						    _dummy_result),
			      digits, _dummy_result);
}

/*! 
 * Compute the integral of the squared entropy error function
 * between a reference solution and the reconstructed
 * polynomial over the domain of cell (iCell,jCell).
 * This result is used in the evaluation of L2 error norm.
 */
template<class SOLN_STATE> 
template<typename Function_Object_Type> inline
double HighOrder2D<SOLN_STATE>::ComputeSolutionEntropyErrorL2(const int &iCell, const int &jCell,
							      Function_Object_Type FuncObj,
							      const int &digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration
  
  return IntegrateOverTheCell(iCell,jCell,
			      square_error_function(FuncObj,
						    wrapped_soln_block_member_function(this,
										       &ClassType::SolutionEntropyAtCoordinates,
										       iCell, jCell,
										       _dummy_result),
						    _dummy_result),
			      digits, _dummy_result);
}

/*! 
 * Compute the integral of the squared error function between the
 * polynomial reconstructions of two high-order variables
 * over the domain of cell (iCell,jCell).
 *
 * \note This subroutine doesn't check that the domains are the same! It's the caller's responsibility.
 */
template<class SOLN_STATE> inline
double HighOrder2D<SOLN_STATE>::ComputeReconstructionsErrorL2(const int &iCell, const int &jCell,
							      const HighOrder2D<Soln_State> & Obj,
							      const unsigned &parameter, const int &digits){

  double _dummy_result(0.0);	 // defined only to provide the return type of the integration
  Vector2D _dummy_Position(0.0); // defined only to provide the type of the position vector where the function is evaluated.

  return IntegrateOverTheCell(iCell,jCell,
			      square_error_function(wrapped_member_function_one_parameter(&Obj,
											  &ClassType::SolutionAtLocation,
											  _dummy_Position,
											  iCell, jCell,
											  parameter,
											  _dummy_result),
						    wrapped_member_function_one_parameter(this,
											  &ClassType::SolutionAtLocation,
											  _dummy_Position,
											  iCell, jCell,
											  parameter,
											  _dummy_result),
						    _dummy_result),
			      digits, _dummy_result);
}

/*! 
 * Calculate the slope limiter for a variety 
 * of limiting functions. 
 */
template<class SOLN_STATE> inline
double HighOrder2D<SOLN_STATE>::CalculateLimiter(double *uQuad, const int &nQuad,
						 const double &u0, const double &u0Min,
						 const double &u0Max, const int &Limiter){
  switch(Limiter) {
  case LIMITER_ONE :
    return ONE;
  case LIMITER_ZERO :
    return ZERO;
  case LIMITER_BARTH_JESPERSEN :
    return Limiter_BarthJespersen(uQuad,u0,u0Min,u0Max,nQuad);
  case LIMITER_VENKATAKRISHNAN :
    return Limiter_Venkatakrishnan(uQuad,u0,u0Min,u0Max,nQuad);
  case LIMITER_VENKATAKRISHNAN_CORRECTED :
    return Limiter_Venkatakrishnan_Modified(uQuad,u0,u0Min,u0Max,nQuad);
  case LIMITER_VANLEER :
    return Limiter_VanLeer(uQuad,u0,u0Min,u0Max,nQuad);
  case LIMITER_VANALBADA :
    return Limiter_VanAlbada(uQuad,u0,u0Min,u0Max,nQuad);
  default:
    return Limiter_BarthJespersen(uQuad,u0,u0Min,u0Max,nQuad);
  };
}

/*! 
 * Output the current object to the 
 * provided output stream.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Output_Object(ostream & out_file) const {

  int i,j, parameter;

  // Output allocation flags
  out_file << _allocated_block << " "
	   << _allocated_cells << " "
	   << _allocated_psinv << " "
	   << _si_calculation  << "\n";

  // Output block indexes
  if (_allocated_block) {
    out_file << Ni << " " << Nj << " " << Ng <<"\n"
	     << OrderOfReconstruction << "\n";
  }

  // Output Taylor derivatives
  if (_allocated_cells){
    for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
      for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) {    
	out_file.setf(ios::skipws,ios::scientific);
	out_file << CellTaylorDeriv(i,j);
	out_file.unsetf(ios::skipws);
	out_file.unsetf(ios::scientific);

	// Output cell inadequate fit flag from last reconstruction
	for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	  out_file << " " << Previous_CellInadequateFitValue(i,j,parameter);
	}

	out_file << endl;
      }
    }
  }
}

/*! 
 * Read the set up of the current object
 * from the provided input stream.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Read_Object(istream & in_file) {

  bool _alloc_block_, _alloc_cells_, _alloc_psinv_;
  int _si_calc_;
  int _Ni_, _Nj_, _Ng_, ReconstructionOrder;
  int i,j, parameter;

  // Read allocation flags
  in_file.setf(ios::skipws);
  in_file >> _alloc_block_
	  >> _alloc_cells_
	  >> _alloc_psinv_
	  >> _si_calc_;

  // Make sure that the execution mode is the same 
  // as at the time when the object was output.
  // Otherwise, the resulting containers will be different.
  if (_si_calc_ != CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS){
    throw runtime_error("HighOrder2D<SOLN_STATE>::Read_Object() ERROR! The object cannot be read due to incompatibilities between the CENO_Execution_Mode class settings and the object settings");
  }

  // check if the object must be allocated
  if (_alloc_block_){
    // Read the block indexes
    in_file >> _Ni_  >> _Nj_  >> _Ng_ 
	    >> ReconstructionOrder;

    // Allocate memory for the object
    allocate(_Ni_-2*_Ng_,
	     _Nj_-2*_Ng_,
	     _Ng_,
	     _alloc_psinv_,
	     ReconstructionOrder);
    
    // check if the cell memory must be allocated
    if (_alloc_cells_){
      // Read the Taylor derivatives
      for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
	for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) { 
	  in_file.setf(ios::skipws);   
	  in_file >> CellTaylorDeriv(i,j);
	  in_file.unsetf(ios::skipws);

	  // Read cell inadequate fit flag for the last reconstruction
	  for (parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	    in_file.setf(ios::skipws);
	    in_file >> Previous_CellInadequateFitValue(i,j,parameter);
	    in_file.unsetf(ios::skipws);
	  }

	  in_file.unsetf(ios::skipws);
	}
      }
      
    } else {
      // Deallocate the cell memory
      deallocate_CellMemory();
    }

  } else {
    // Deallocate the current object
    deallocate();
  }

  in_file.unsetf(ios::skipws);
}

/*! 
 * Output the reconstruction type of each cell
 * to the provided output stream.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::outputReconstructionTypeMap(ostream & out_file) const{

  int i,j;

  if (ReconstructionTypeMap != NULL){
    out_file << endl;
    for (j = Nj-1; j>=0; --j){
      for (i=0; i<Ni; ++i){
	out_file << getCellReconstructionType(i,j);
      }
      out_file << endl;
    }
  } else {
    out_file << "NULL\n";
  }
}

/*!
 * Broadcast high-order object to all            
 * processors involved in the calculation from the      
 * primary processor using the MPI broadcast routine.
 * 
 * \param Block_Geometry the grid block associated to the high-order object
 *
 * \note In the current implementation, only the container
 *       of TaylorDerivatives needs to be broadcast.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Broadcast_HighOrder_Data(GeometryType & Block_Geometry){

#ifdef _MPI_VERSION
  
  int i, j, buffer_size, TD_Bcast, td, var, counter;
  double *buffer;

  // Obs: It is assumed that InitializeBasicVariable() routine
  //      has been already run on each CPU.

  /* Calculate how many derivatives must be broadcast, including the limiter value for each parameter */
  TD_Bcast = NumberOfVariables() * (NumberOfTaylorDerivatives() + 1 );

  /* Calculate the size of the buffer.*/
  buffer_size = TD_Bcast*Ni*Nj;
  
  buffer = new double[buffer_size];

  // Load the buffer on the primary CPU
  if (CFFC_Primary_MPI_Processor()) {
    counter = 0;
    // Pack all the information 
    for (j  = 0 ; j < Nj ; ++j ) { //< for each jCell
      for ( i = 0 ; i < Ni ; ++i ) {  //< for each iCell
	for (var = 1; var <= NumberOfVariables(); ++var){ //< for each state parameter

	  // Load the derivatives of the current solution state parameter
	  for (td = 0; td <= CellTaylorDeriv(i,j).LastElem(); ++td, ++counter){ //< for each Taylor derivative of (iCell,jCell)
	    buffer[counter] = CellTaylorDerivState(i,j,td)[var];
	  } /* endfor(td) */

	  // Load the limiter of the current solution state parameter
	  buffer[counter] = CellTaylorDeriv(i,j).Limiter(var);
	  ++counter;
	}/* endfor(var) */
      } /* endfor(i) */
    } /* endfor(j) */

  }/* endif(CFFC_Primary_MPI_Processor()) */


  // Broadcast buffer
  MPI::COMM_WORLD.Bcast(buffer, buffer_size, MPI::DOUBLE, 0);

  // Unload the buffer on the receiver CPUs and set the variables
  if (!CFFC_Primary_MPI_Processor()) {
    counter = 0;
    // Unpack all the information 
    for (j  = 0 ; j < Nj; ++j ) { //< for each jCell
      for ( i = 0; i < Ni; ++i ) {  //< for each iCell
  	for (var = 1; var <= NumberOfVariables(); ++var){ //< for each state parameter

	  // Unload the derivatives of the current solution state parameter
	  for (td = 0; td <= CellTaylorDeriv(i,j).LastElem(); ++td, ++counter){ //< for each Taylor derivative of (iCell,jCell)
	    CellTaylorDerivState(i,j,td)[var] = buffer[counter];
	  } /* endfor(td) */
	  
	  // Unload the limiter of the current solution state parameter
	  CellTaylorDeriv(i,j).Limiter(var) = buffer[counter];
	  ++counter;

	}/* endfor(var) */
      } /* endfor(i) */
    } /* endfor(j) */
  }/* endif(!CFFC_Primary_MPI_Processor()) */


  // Deallocate memory
  delete []buffer;
  buffer = NULL;

#endif

}

#ifdef _MPI_VERSION
/*!
 * Broadcast high-order object to all processors 
 * associated with the specified communicator from the  
 * specified processor using the MPI broadcast routine.
 *
 * \param Communicator a particular MPI communicator
 * \param Source_CPU the CPU used as source for the broadcast
 * \param Block_Geometry the grid block to which the
 *                   high-order object is associated. 
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Broadcast_HighOrder_Data(MPI::Intracomm &Communicator, 
						       const int &Source_CPU,
						       GeometryType & Block_Geometry){

  int Source_Rank = 0; 
  int i, j, buffer_size, TD_Bcast, td, var, counter;
  double *buffer;
  
  // Obs: It is assumed that InitializeBasicVariable() routine
  //      has been already run on each involved CPU.

  /* Calculate how many derivatives must be broadcast, including the limiter value for each parameter */
  TD_Bcast = NumberOfVariables() * (NumberOfTaylorDerivatives() + 1 );

  /* Calculate the size of the buffer.*/
  buffer_size = TD_Bcast*Ni*Nj;
  
  buffer = new double[buffer_size];

  // Load the buffer on the source CPU
  if (CFFC_MPI::This_Processor_Number == Source_CPU) {
    counter = 0;
    // Pack all the information 
    for (j  = 0 ; j < Nj ; ++j ) { //< for each jCell
      for ( i = 0 ; i < Ni ; ++i ) {  //< for each iCell
	for (var = 1; var <= NumberOfVariables(); ++var){ //< for each state parameter

	  // Load the derivatives of the current solution state parameter
	  for (td = 0; td <= CellTaylorDeriv(i,j).LastElem(); ++td, ++counter){ //< for each Taylor derivative of (iCell,jCell)
	    buffer[counter] = CellTaylorDerivState(i,j,td)[var];
	  } /* endfor(td) */

	  // Load the limiter of the current solution state parameter
	  buffer[counter] = CellTaylorDeriv(i,j).Limiter(var);
	  ++counter;
	}/* endfor(var) */
      } /* endfor(i) */
    } /* endfor(j) */

  }/* endif(CFFC_MPI::This_Processor_Number == Source_CPU) */


  // Broadcast buffer
  Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_Rank);

  // Unload the buffer on the receiver CPUs and set the variables
  if (CFFC_MPI::This_Processor_Number != Source_CPU) {
    counter = 0;
    // Unpack all the information 
    for (j  = 0 ; j < Nj; ++j ) { //< for each jCell
      for ( i = 0; i < Ni; ++i ) {  //< for each iCell
  	for (var = 1; var <= NumberOfVariables(); ++var){ //< for each state parameter

	  // Unload the derivatives of the current solution state parameter
	  for (td = 0; td <= CellTaylorDeriv(i,j).LastElem(); ++td, ++counter){ //< for each Taylor derivative of (iCell,jCell)
	    CellTaylorDerivState(i,j,td)[var] = buffer[counter];
	  } /* endfor(td) */
	  
	  // Unload the limiter of the current solution state parameter
	  CellTaylorDeriv(i,j).Limiter(var) = buffer[counter];
	  ++counter;

	}/* endfor(var) */
      } /* endfor(i) */
    } /* endfor(j) */
  }/* endif(CFFC_MPI::This_Processor_Number != Source_CPU) */


  // Deallocate memory
  delete []buffer;
  buffer = NULL;  
  
}
#endif

/*!
 * Compute the AMR criteria for the block based on the 
 * minimum smoothness indicator value encountered over the
 * block cells and all solution variables.
 * \todo Add more comments here!
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
double HighOrder2D<SOLN_STATE>::AMR_Criteria_Based_On_Minimum_Smoothness_Indicator(Soln_Block_Type &SolnBlk){

  if (CENO_Tolerances::Fit_Tolerance <= ZERO){
    throw runtime_error("HighOrder2D<SOLN_STATE>::AMR_Criteria_Based_On_Minimum_Smoothness_Indicator() ERROR! Negative/zero CENO tolerance is not allowed for refinement!");
  }

  int i, j, parameter;

  double SI_Min;		//< minimum smoothness indicator value

  // Reconstruct the block solution (i.e. high-order, data analysis and monotonicity enforcement).
  ComputeHighOrderSolutionReconstruction(SolnBlk,
					 CENO_Execution_Mode::Limiter);

  /* Initialize the minimum smoothness indicator value for the block. */
  SI_Min = CellSmoothnessIndicatorValue(ICl,JCl,1);

  /* Calculate the minimum smoothness indicator value for 
     all block interior cells and solution variables. */
  for ( j  = JCl ; j <= JCu ; ++j ) {
    for ( i = ICl ; i <= ICu ; ++i ) {
      for(parameter = 1; parameter <= NumberOfVariables(); ++parameter){
	SI_Min = min(SI_Min, CellSmoothnessIndicatorValue(i,j,parameter));
      } /* endfor (parameter) */
    } /* endfor (i) */
  } /* endfor (j) */

  // Evaluate the block refinement criteria based on the minimum encountered SI.
  return ( exp( -max(ZERO, SI_Min) / (CENO_Tolerances::AMR_Smoothness_Units * CENO_Tolerances::Fit_Tolerance )) );
  
}

/*!
 * Build the map of reconstruction types,
 * which means to identify the type of reconstruction
 * of each cell in a constrained block.
 * This routine assumes that the block type has been 
 * already identified and the properties of the 
 * high-order block boundaries have been already setup.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::BuildReconstructionTypeMap(void){

  int i,j;

  // Allocate memory for the ReconstructionTypeMap variable
  allocate_ReconstructionTypeMap();

  if (IsConstrainedReconstructionRequired()){

    // Set the NO reconstruction
    for (j = 0; j < Nj; ++j){
      for (i = 0; i < ICl-Nghost_HO; ++i){
	ReconstructionTypeMap[i][j] = 'n';
      }

      for (i = ICu+Nghost_HO+1; i < Ni; ++i){
	ReconstructionTypeMap[i][j] = 'n';
      }
    }

    for (i = ICl-Nghost_HO; i <= ICu+Nghost_HO; ++i){
      for (j = 0; j < JCl-Nghost_HO; ++j){
	ReconstructionTypeMap[i][j] = 'n';
      }

      for (j = JCu+Nghost_HO+1; j < Nj; ++j){
	ReconstructionTypeMap[i][j] = 'n';
      }
    }

    // Set regular reconstruction where central reconstruction normally occurs in the absence of constrained boundaries
    for (i = ICl-Nghost_HO; i <= ICu+Nghost_HO; ++i){
      for (j = JCl-Nghost_HO; j <= JCu+Nghost_HO; ++j){
	ReconstructionTypeMap[i][j] = 'r';
      }
    }

    // ==== North Boundary ====
    if (NorthBnd.IsReconstructionConstrained()){
      for (i = ICl; i <= ICu; ++i){
	// Set the NO reconstruction cells
	for (j = JCu+1; j<=JCu+Nghost_HO; ++j){
	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[i][JCu] = 'c';
      
	// Set the unconstrained cells with modified stencil
	for (j = JCu-1; j>JCu-rings; --j){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }
    }

    // ==== South Boundary ====
    if (SouthBnd.IsReconstructionConstrained()){
      for (i = ICl; i <= ICu; ++i){
	// Set the NO reconstruction cells
	for (j = JCl-Nghost_HO; j<JCl; ++j){
	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[i][JCl] = 'c';
      
	// Set the unconstrained cells with modified stencil
	for (j = JCl+1; j<JCl+rings; ++j){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }
    }

    // ==== West Boundary ====
    if (WestBnd.IsReconstructionConstrained()){
      for (j = JCl; j<=JCu; ++j){
	// Set the NO reconstruction cells
	for (i=ICl-Nghost_HO; i<ICl; ++i){
	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells 
	ReconstructionTypeMap[ICl][j] = 'c';
      
	// Set the unconstrained cells with modified stencil
	for (i=ICl+1; i<ICl+rings; ++i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      // Check the North and South boundaries
      if (NorthBnd.IsReconstructionConstrained()){
	// Fix the cells
	for (i=ICl+1; i<ICl+rings; ++i){
	  ReconstructionTypeMap[i][JCu] = 'c';
	}
      }

      if (SouthBnd.IsReconstructionConstrained()){
	// Fix the cells
	for (i=ICl+1; i<ICl+rings; ++i){
	  ReconstructionTypeMap[i][JCl] = 'c';
	}
      }
    }

    // ==== East Boundary ====
    if (EastBnd.IsReconstructionConstrained()){
      for (j = JCl; j<=JCu; ++j){
	// Set the NO reconstruction cells
	for (i=ICu+1; i<=ICu+Nghost_HO; ++i){
	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells 
	ReconstructionTypeMap[ICu][j] = 'c';
      
	// Set the unconstrained cells with modified stencil
	for (i=ICu-1; i>ICu-rings; --i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      // Check the North and South boundaries
      if (NorthBnd.IsReconstructionConstrained()){
	// Fix the cells
	for (i=ICu-1; i>ICu-rings; --i){
	  ReconstructionTypeMap[i][JCu] = 'c';
	}
      }

      if (SouthBnd.IsReconstructionConstrained()){
	// Fix the cells
	for (i=ICu-1; i>ICu-rings; --i){
	  ReconstructionTypeMap[i][JCl] = 'c';
	}
      }
    }

    // ==== Corner North-West ====
    if (N_WestBnd.IsReconstructionConstrained()){
      for (j = JCu+1; j<= JCu+Nghost_HO; ++j){
	// Set the NO reconstruction cells
	for (i=ICl-Nghost_HO; i<=ICl-1; ++i){
 	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[ICl][j] = 'c';

	// Set the unconstrained cells with modified stencil
	for (i = ICl+1; i<ICl+rings; ++i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      if (!WestBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICl-Nghost_HO; i<ICl+rings; ++i){
	  for (j = JCu; j > JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else if (N_WestBnd.IsReconstructionStencilAffected()){

      if (W_NorthBnd.IsReconstructionConstrained()){

	for (i = ICl-Nghost_HO; i<ICl; ++i){
	  // Set the NO reconstruction cells
	  for (j = JCu+1; j<=JCu+Nghost_HO; ++j){
	    ReconstructionTypeMap[i][j] = 'n';
	  }

	  // Set the constrained cells
	  ReconstructionTypeMap[i][JCu] = 'c';

	  // Set the unconstrained cells with modified stencil
	  for (j = JCu-1; j>JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}

	if (!NorthBnd.IsReconstructionStencilAffected()){
	  // Set the cells with modified stencil
	  for (i = ICl; i<ICl+rings; ++i){
	    for (j = JCu+Nghost_HO; j > JCu-rings; --j){
	      ReconstructionTypeMap[i][j] = 'm';
	    }
	  }
	} 

      } else if (W_NorthBnd.IsReconstructionStencilAffected()) {
	// Both extensions affect the stencil. 
	// Set the corner ghost cells to NO reconstruction
	for (j = JCu+1; j<= JCu+Nghost_HO; ++j){
	  for (i=ICl-Nghost_HO; i<=ICl-1; ++i){
	    ReconstructionTypeMap[i][j] = 'n';
	  }
	}

      } else {
	// W_NorthBnd is transparent
	
	// Set the cells with modified stencil
	for (i = ICl-1; i>=ICl-rings; --i){
	  for (j = JCu+Nghost_HO; j > JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else {
      // N_WestBnd is transparent

      if (W_NorthBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICl-Nghost_HO; i<ICl+rings; ++i){
	  for (j = JCu+1; j<=JCu+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } // endif (N_WestBnd)


    // ==== Corner South-West ====
    if (S_WestBnd.IsReconstructionConstrained()){
      for (j = JCl-1; j>=JCl-Nghost_HO; --j){
	// Set the NO reconstruction cells
	for (i=ICl-Nghost_HO; i<ICl; ++i){
 	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[ICl][j] = 'c';

	// Set the unconstrained cells with modified stencil
	for (i = ICl+1; i<ICl+rings; ++i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      if (!WestBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICl-Nghost_HO; i<ICl+rings; ++i){
	  for (j = JCl; j < JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else if (S_WestBnd.IsReconstructionStencilAffected()){

      if (W_SouthBnd.IsReconstructionConstrained()){
	
	for (i = ICl-Nghost_HO; i<ICl; ++i){
	  // Set the NO reconstruction cells
	  for (j = JCl-1; j>=JCl-Nghost_HO; --j){
	    ReconstructionTypeMap[i][j] = 'n';
	  }

	  // Set the constrained cells
	  ReconstructionTypeMap[i][JCl] = 'c';

	  // Set the unconstrained cells with modified stencil
	  for (j = JCl+1; j<JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}

	if (!SouthBnd.IsReconstructionStencilAffected()){
	  // Set the cells with modified stencil
	  for (i = ICl; i<ICl+rings; ++i){
	    for (j = JCl-Nghost_HO; j < JCl+rings; ++j){
	      ReconstructionTypeMap[i][j] = 'm';
	    }
	  }
	}

      } else if (W_SouthBnd.IsReconstructionStencilAffected()) {
	// Both extensions affect the stencil. 
	// Set the corner ghost cells to NO reconstruction
	for (j = JCl-Nghost_HO; j<= JCl-1; ++j){
	  for (i=ICl-Nghost_HO; i<=ICl-1; ++i){
	    ReconstructionTypeMap[i][j] = 'n';
	  }
	}

      } else {
	// W_SouthBnd is transparent
	
	// Set the cells with modified stencil
	for (i = ICl-1; i>=ICl-rings; --i){
	  for (j = JCl-Nghost_HO; j < JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else {
      // S_WestBnd is transparent

      if (W_SouthBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICl-Nghost_HO; i<ICl+rings; ++i){
	  for (j = JCl-1; j>=JCl-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }
      
    } // endif (S_WestBnd)


    // ==== Corner South-East ====
    if (S_EastBnd.IsReconstructionConstrained()){

      for (j = JCl-1; j>=JCl-Nghost_HO; --j){
	// Set the NO reconstruction cells
	for (i=ICu+1; i<=ICu+Nghost_HO; ++i){
 	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[ICu][j] = 'c';

	// Set the unconstrained cells with modified stencil
	for (i = ICu-1; i>ICu-rings; --i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      if (!EastBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICu+Nghost_HO; i>ICu-rings; --i){
	  for (j = JCl; j < JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else if (S_EastBnd.IsReconstructionStencilAffected()){
      
      if (E_SouthBnd.IsReconstructionConstrained()){
	
	for (i = ICu+1; i<=ICu+Nghost_HO; ++i){
	  // Set the NO reconstruction cells
	  for (j = JCl-1; j>=JCl-Nghost_HO; --j){
	    ReconstructionTypeMap[i][j] = 'n';
	  }

	  // Set the constrained cells
	  ReconstructionTypeMap[i][JCl] = 'c';

	  // Set the unconstrained cells with modified stencil
	  for (j = JCl+1; j<JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}

	if (!SouthBnd.IsReconstructionStencilAffected()){
	  // Set the cells with modified stencil
	  for (i = ICu; i>ICu-rings; --i){
	    for (j = JCl-Nghost_HO; j < JCl+rings; ++j){
	      ReconstructionTypeMap[i][j] = 'm';
	    }
	  }
	}

      } else if (E_SouthBnd.IsReconstructionStencilAffected()) {
	// Both extensions affect the stencil. 
	// Set the corner ghost cells to NO reconstruction
	for (j = JCl-Nghost_HO; j<= JCl-1; ++j){
	  for (i=ICu+1; i<=ICu+Nghost_HO; ++i){
	    ReconstructionTypeMap[i][j] = 'n';
	  }
	}

      } else {
	// E_SouthBnd is transparent
	
	// Set the cells with modified stencil
	for (i = ICu+1; i<=ICu+rings; ++i){
	  for (j = JCl-Nghost_HO; j < JCl+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else {
      // S_EastBnd is transparent

      if (E_SouthBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICu+Nghost_HO; i>ICu-rings; --i){
	  for (j = JCl-1; j>=JCl-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } // endif (S_EastBnd)


    // ==== Corner North-East ====
    if (N_EastBnd.IsReconstructionConstrained()){

      for (j = JCu+1; j<= JCu+Nghost_HO; ++j){
	// Set the NO reconstruction cells
	for (i=ICu+1; i<=ICu+Nghost_HO; ++i){
 	  ReconstructionTypeMap[i][j] = 'n';
	}

	// Set the constrained cells
	ReconstructionTypeMap[ICu][j] = 'c';

	// Set the unconstrained cells with modified stencil
	for (i = ICu-1; i>ICu-rings; --i){
	  ReconstructionTypeMap[i][j] = 'm';
	}
      }

      if (!EastBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICu+Nghost_HO; i>ICu-rings; --i){
	  for (j = JCu; j > JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else if (N_EastBnd.IsReconstructionStencilAffected()){

      if (E_NorthBnd.IsReconstructionConstrained()){

	for (i = ICu+1; i<=ICu+Nghost_HO; ++i){
	  // Set the NO reconstruction cells
	  for (j = JCu+1; j<=JCu+Nghost_HO; ++j){
	    ReconstructionTypeMap[i][j] = 'n';
	  }

	  // Set the constrained cells
	  ReconstructionTypeMap[i][JCu] = 'c';

	  // Set the unconstrained cells with modified stencil
	  for (j = JCu-1; j>JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}

	if (!NorthBnd.IsReconstructionStencilAffected()){
	  // Set the cells with modified stencil
	  for (i = ICu; i>ICu-rings; --i){
	    for (j = JCu+Nghost_HO; j > JCu-rings; --j){
	      ReconstructionTypeMap[i][j] = 'm';
	    }
	  }
	}

      } else if (E_NorthBnd.IsReconstructionStencilAffected()) {
	// Both extensions affect the stencil. 
	// Set the corner ghost cells to NO reconstruction
	for (j = JCu+1; j<= JCu+Nghost_HO; ++j){
	  for (i=ICu+1; i<=ICu+Nghost_HO; ++i){
	    ReconstructionTypeMap[i][j] = 'n';
	  }
	}

      } else {
	// E_NorthBnd is transparent
	
	// Set the cells with modified stencil
	for (i = ICu+1; i<=ICu+rings; ++i){
	  for (j = JCu+Nghost_HO; j > JCu-rings; --j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } else {
      // N_EastBnd is transparent

      if (E_NorthBnd.IsReconstructionStencilAffected()){
	// Set the cells with modified stencil
	for (i = ICu+Nghost_HO; i>ICu-rings; --i){
	  for (j = JCu+1; j<=JCu+rings; ++j){
	    ReconstructionTypeMap[i][j] = 'm';
	  }
	}
      }

    } // endif (N_EastBnd)

  } else {
    // The block doesn't have ReconstructionTypeMap variable allocated
    return;
  }
}

// Friend functions

//! operator<<
template<class SOLN_STATE> inline
ostream & operator<< (ostream & os, const HighOrder2D<SOLN_STATE> & Obj){
  Obj.Output_Object(os);
  return os;
}

//! operator>>
template<class SOLN_STATE> inline
istream & operator>> (istream & os, HighOrder2D<SOLN_STATE> & Obj){
  Obj.Read_Object(os);
  return os;
}


        /*:::::::::::::::::::*/                                             
        /*  Specializations  */                                             
        /*:::::::::::::::::::*/     

//! Number of variables for type 'double'
template<> inline
int HighOrder2D<double>::NumberOfVariables(void) const {
  return 1;
}

//! Get the value of Taylor derivative of cell (ii,jj) for the (p1,p2) powers and the specified Variable for type 'double'
template<> inline
const double & HighOrder2D<double>::CellTaylorDerivValue(const int & ii, const int & jj,
							 const int & p1, const int & p2, const int & Variable) const {
  return TD[ii][jj](p1,p2);
}

//! Get the value of Taylor derivative of cell (ii,jj) for the (p1,p2) powers and the specified Variable for type 'double'
template<> inline
double & HighOrder2D<double>::CellTaylorDerivValue(const int & ii, const int & jj,
						   const int & p1, const int & p2, const int & Variable) {
  return TD[ii][jj](p1,p2);
}

template <> inline
double HighOrder2D<double>::SolutionAtCoordinates(const int & ii, const int & jj, 
						  const double & X_Coord, const double & Y_Coord,
						  const unsigned & parameter) const {
  return TD[ii][jj].ComputeSolutionFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
}


/* ---------------------------------------------------------------------------------------------- 
 * =============== INCLUDE THE IMPLEMENTATION OF HIGH-ORDER 2D RECONSTRUCTIONS ==================
 * ---------------------------------------------------------------------------------------------*/
#include "HighOrder2D_Reconstructions.h"

#endif // _HIGHORDER_2D_INCLUDED
