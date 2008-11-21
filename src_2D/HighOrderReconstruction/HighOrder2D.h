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

  //! @name Indexes of cells between which the smoothness indicator is computed with central stencil
  //@{
  const int & StartIdir_SI(void) const { return StartI_SI; }
  const int & EndIdir_SI(void) const { return EndI_SI; }
  const int & StartJdir_SI(void) const { return StartJ_SI; }
  const int & EndJdir_SI(void) const { return EndJ_SI; }
  //@}

  //! @name Indexes of cells between which the limited piecewise linear (LPWL) reconstruction is potentially applied.
  //@{
  const int & StartIdir_LPWL(void) const { return StartI_LPWL; }
  const int & EndIdir_LPWL(void) const { return EndI_LPWL; }
  const int & StartJdir_LPWL(void) const { return StartJ_LPWL; }
  const int & EndJdir_LPWL(void) const { return EndJ_LPWL; }
  //@}

  //! @name Indexes of cells between which constrained reconstruction is calculated if West BC is imposed by constraints.
  //@{
  const int & StartIdir_ConstrWest(void) const { return StartI_ConstrWest; }
  const int & EndIdir_ConstrWest(void) const { return EndI_ConstrWest; }
  const int & StartJdir_ConstrWest(void) const { return StartJ_ConstrWest; }
  const int & EndJdir_ConstrWest(void) const { return EndJ_ConstrWest; }
  //@}

  //! @name Indexes of cells between which constrained reconstruction is calculated if East BC is imposed by constraints.
  //@{
  const int & StartIdir_ConstrEast(void) const { return StartI_ConstrEast; }
  const int & EndIdir_ConstrEast(void) const { return EndI_ConstrEast; }
  const int & StartJdir_ConstrEast(void) const { return StartJ_ConstrEast; }
  const int & EndJdir_ConstrEast(void) const { return EndJ_ConstrEast; }
  //@}

  //! @name Indexes of cells between which constrained reconstruction is calculated if North BC is imposed by constraints.
  //@{
  const int & StartIdir_ConstrNorth(void) const { return StartI_ConstrNorth; }
  const int & EndIdir_ConstrNorth(void) const { return EndI_ConstrNorth; }
  const int & StartJdir_ConstrNorth(void) const { return StartJ_ConstrNorth; }
  const int & EndJdir_ConstrNorth(void) const { return EndJ_ConstrNorth; }
  //@}

  //! @name Indexes of cells between which constrained reconstruction is calculated if South BC is imposed by constraints.
  //@{
  const int & StartIdir_ConstrSouth(void) const { return StartI_ConstrSouth; }
  const int & EndIdir_ConstrSouth(void) const { return EndI_ConstrSouth; }
  const int & StartJdir_ConstrSouth(void) const { return StartJ_ConstrSouth; }
  const int & EndJdir_ConstrSouth(void) const { return EndJ_ConstrSouth; }
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
  const bool & IsWestConstrainedReconstructionRequired(void) const { return _constrained_WEST_reconstruction; }
  //! Return true if constrained reconstruction is used in the proximity of East block boundary, otherwise return false
  const bool & IsEastConstrainedReconstructionRequired(void) const { return _constrained_EAST_reconstruction; }
  //! Return true if constrained reconstruction is used in the proximity of North block boundary, otherwise return false
  const bool & IsNorthConstrainedReconstructionRequired(void) const { return _constrained_NORTH_reconstruction; }
  //! Return true if constrained reconstruction is used in the proximity of South block boundary, otherwise return false
  const bool & IsSouthConstrainedReconstructionRequired(void) const { return _constrained_SOUTH_reconstruction; }
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
  Soln_State XGradientAtCoordinates(const int & ii, const int & jj,
				    const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeXGradientFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }
  //! Evaluate the X-gradient of the interpolant at a given position vector for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State XGradientStateAtLocation(const int & ii, const int & jj, const Vector2D &CalculationPoint) const {
    return XGradientAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y);
  }

  //! Evaluate the Y-gradient of the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State YGradientAtCoordinates(const int & ii, const int & jj,
				    const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeYGradientFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }
  //! Evaluate the Y-gradient of the interpolant at a given position vector for all solution variables,
  //  using the reconstruction of cell (ii,jj).
  Soln_State YGradientStateAtLocation(const int & ii, const int & jj, const Vector2D &CalculationPoint) const {
    return YGradientAtCoordinates(ii,jj,CalculationPoint.x,CalculationPoint.y);
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
  void ComputeCellReconstructionPseudoInverseNearConstrainedBoundaries(const int& BOUNDARY,
								       const int& iCell, const int& jCell);

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
  /*! @brief Set the range of cells with straight edges in which reconstruction is performed.  */
  void SetRangeOfQuadCellsWithoutConstrainedReconstruction(void);
  /*! @brief Set the range of cells in which constrained reconstruction is performed.  */
  void SetRangeOfQuadCellsWithConstrainedReconstruction(void);
  /*! @brief Set the central stencil of cells used for reconstruction.  */
  void SetReconstructionStencil(const int &iCell, const int &jCell,
				IndexType & i_index, IndexType & j_index) const;
  void SetSpecialReconstructionStencil(const int &iCell, const int &jCell,
				       IndexType & i_index, IndexType & j_index) const;
  void SetConstrainedReconstructionStencil(const int &iCell, const int &jCell,
					   IndexType & i_index, IndexType & j_index) const;  
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
  int StartI_SI,               //!< Index of the first cell in i-direction in which SI is computed with centered stencil
    EndI_SI,                   //!< Index of the last cell in i-direction in which SI is computed with centered stencil
    StartJ_SI,                 //!< Index of the first cell in j-direction in which SI is computed with centered stencil
    EndJ_SI;                   //!< Index of the last cell in j-direction in which SI is computed with centered stencil
  int StartI_LPWL,             //!< Index of the first cell in i-direction which is checked for non-smooth interpolant
    EndI_LPWL,                 //!< Index of the last cell in i-direction which is checked for non-smooth interpolant
    StartJ_LPWL,               //!< Index of the first cell in j-direction which is checked for non-smooth interpolant
    EndJ_LPWL;                 //!< Index of the last cell in j-direction which is checked for non-smooth interpolant
  int StartI_ConstrWest,       //!< Index of the first cell in i-direction which is influenced by a constrained West boundary 
    EndI_ConstrWest,           //!< Index of the last cell in i-direction which is influenced by a constrained West boundary
    StartJ_ConstrWest,         //!< Index of the first cell in j-direction which is influenced by a constrained West boundary
    EndJ_ConstrWest;           //!< Index of the last cell in j-direction which is influenced by a constrained West boundary
  int StartI_ConstrEast,       //!< Index of the first cell in i-direction which is influenced by a constrained East boundary 
    EndI_ConstrEast,           //!< Index of the last cell in i-direction which is influenced by a constrained East boundary
    StartJ_ConstrEast,         //!< Index of the first cell in j-direction which is influenced by a constrained East boundary
    EndJ_ConstrEast;           //!< Index of the last cell in j-direction which is influenced by a constrained East boundary
  int StartI_ConstrNorth,      //!< Index of the first cell in i-direction which is influenced by a constrained North boundary 
    EndI_ConstrNorth,          //!< Index of the last cell in i-direction which is influenced by a constrained North boundary
    StartJ_ConstrNorth,        //!< Index of the first cell in j-direction which is influenced by a constrained North boundary
    EndJ_ConstrNorth;          //!< Index of the last cell in j-direction which is influenced by a constrained North boundary
  int StartI_ConstrSouth,      //!< Index of the first cell in i-direction which is influenced by a constrained South boundary 
    EndI_ConstrSouth,          //!< Index of the last cell in i-direction which is influenced by a constrained South boundary
    StartJ_ConstrSouth,        //!< Index of the first cell in j-direction which is influenced by a constrained South boundary
    EndJ_ConstrSouth;          //!< Index of the last cell in j-direction which is influenced by a constrained South boundary
  int StartJ_SI_ConstrWest,    //!< Index of the first cell in j-direction influenced by a constrained West boundary for SI calc.
    EndJ_SI_ConstrWest;        //!< Index of the last cell in j-direction influenced by a constrained West boundary for SI calc.
  int StartJ_SI_ConstrEast,    //!< Index of the first cell in j-direction influenced by a constrained East boundary for SI calc.
    EndJ_SI_ConstrEast;        //!< Index of the last cell in j-direction influenced by a constrained East boundary for SI calc.
  int StartI_SI_ConstrNorth,   //!< Index of the first cell in i-direction influenced by a constrained North boundary for SI calc.
    EndI_SI_ConstrNorth;       //!< Index of the last cell in i-direction influenced by a constrained North boundary for SI calc.
  int StartI_SI_ConstrSouth,   //!< Index of the first cell in i-direction influenced by a constrained South boundary for SI calc.
    EndI_SI_ConstrSouth;       //!< Index of the last cell in i-direction influenced by a constrained South boundary for SI calc.
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
  bool _constrained_WEST_reconstruction; //!< Flag for the WEST boundary to show if constrained rec. is required.
  bool _constrained_EAST_reconstruction; //!< Flag for the EAST boundary to show if constrained rec. is required.
  bool _constrained_NORTH_reconstruction; //!< Flag for the NORTH boundary to show if constrained rec. is required.
  bool _constrained_SOUTH_reconstruction; //!< Flag for the SOUTH boundary to show if constrained rec. is required.
  //@}

  //! Get the number of variables in the solution state
  int NumberOfVariables(void) const {return Soln_State::NumVar(); }

  //! Allocate memory at the cell level based on the order of reconstruction
  void allocate_CellMemory(const int &ReconstructionOrder, const bool &_pseudo_inverse_allocation_);

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

  /*! @brief Set the central stencil of cells used to compute the smoothness indicator.  */
  void SetSmoothnessIndicatorStencil(const int &iCell, const int &jCell);
  //! @brief Flag (iCell,jCell) cell with non-smooth high-order interpolants.
  void FlagCellReconstructionsAsNonSmooth(const int &iCell, const int &jCell);
  /*! @brief Set the stencil of cells used to compute the smoothness indicator for a cell affected by constrained boundaries.*/
  void SetSmoothnessIndicatorStencilForConstrainedReconstruction(const int &iCell, const int &jCell);
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
  StartI_SI(0), EndI_SI(0), StartJ_SI(0), EndJ_SI(0),
  StartI_LPWL(0), EndI_LPWL(0), StartJ_LPWL(0), EndJ_LPWL(0),
  StartI_ConstrWest(0), EndI_ConstrWest(0), StartJ_ConstrWest(0), EndJ_ConstrWest(0),
  StartI_ConstrNorth(0), EndI_ConstrNorth(0), StartJ_ConstrNorth(0), EndJ_ConstrNorth(0),
  StartI_ConstrSouth(0), EndI_ConstrSouth(0), StartJ_ConstrSouth(0), EndJ_ConstrSouth(0),
  StartI_ConstrEast(0), EndI_ConstrEast(0), StartJ_ConstrEast(0), EndJ_ConstrEast(0),
  StartJ_SI_ConstrWest(0),  EndJ_SI_ConstrWest(0),
  StartI_SI_ConstrNorth(0), EndI_SI_ConstrNorth(0),
  StartI_SI_ConstrSouth(0), EndI_SI_ConstrSouth(0),
  StartJ_SI_ConstrEast(0),  EndJ_SI_ConstrEast(0),
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  _freeze_limiter(false),
  TD(NULL), SI(NULL), LimitedCell(NULL), PreviousLimitedCell(NULL),
  rings(0), rings_SI(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL),
  Geom(NULL), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS),
  _constrained_block_reconstruction(false),
  _constrained_WEST_reconstruction(false) , _constrained_EAST_reconstruction(false),
  _constrained_NORTH_reconstruction(false), _constrained_SOUTH_reconstruction(false),
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
  StartI_SI(0), EndI_SI(0), StartJ_SI(0), EndJ_SI(0),
  StartI_LPWL(0), EndI_LPWL(0), StartJ_LPWL(0), EndJ_LPWL(0),
  StartI_ConstrWest(0), EndI_ConstrWest(0), StartJ_ConstrWest(0), EndJ_ConstrWest(0),
  StartI_ConstrNorth(0), EndI_ConstrNorth(0), StartJ_ConstrNorth(0), EndJ_ConstrNorth(0),
  StartI_ConstrSouth(0), EndI_ConstrSouth(0), StartJ_ConstrSouth(0), EndJ_ConstrSouth(0),
  StartI_ConstrEast(0), EndI_ConstrEast(0), StartJ_ConstrEast(0), EndJ_ConstrEast(0),
  StartJ_SI_ConstrWest(0),  EndJ_SI_ConstrWest(0),
  StartI_SI_ConstrNorth(0), EndI_SI_ConstrNorth(0),
  StartI_SI_ConstrSouth(0), EndI_SI_ConstrSouth(0),
  StartJ_SI_ConstrEast(0),  EndJ_SI_ConstrEast(0),
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  _freeze_limiter(false),
  TD(NULL), SI(NULL), LimitedCell(NULL), PreviousLimitedCell(NULL),
  rings(0), rings_SI(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL), 
  Geom(&Block), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS),
  _constrained_block_reconstruction(false),
  _constrained_WEST_reconstruction(false) , _constrained_EAST_reconstruction(false),
  _constrained_NORTH_reconstruction(false), _constrained_SOUTH_reconstruction(false),
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

  // Determine the range of cells in which NO constrained reconstruction is performed
  SetRangeOfQuadCellsWithoutConstrainedReconstruction();

  // Determine the range of cells in which !constrained! reconstruction is performed
  SetRangeOfQuadCellsWithConstrainedReconstruction();
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

    // set range of cells with central stencil smoothness indicator calculation
    StartI_SI = rhs.StartI_SI;
    EndI_SI = rhs.EndI_SI;
    StartJ_SI = rhs.StartJ_SI;
    EndJ_SI = rhs.EndJ_SI;

    // set range of cells checked for non-smooth interpolant
    StartI_LPWL = rhs.StartI_LPWL;
    EndI_LPWL = rhs.EndI_LPWL;
    StartJ_LPWL = rhs.StartJ_LPWL;
    EndJ_LPWL = rhs.EndJ_LPWL;

    // set range of cells which have constrained reconstruction if the West boundary condition is enforced with constraints
    StartI_ConstrWest = rhs.StartI_ConstrWest;
    EndI_ConstrWest = rhs.EndI_ConstrWest;
    StartJ_ConstrWest = rhs.StartJ_ConstrWest;
    EndJ_ConstrWest = rhs.EndJ_ConstrWest;

    // set range of cells which have constrained reconstruction if the East boundary condition is enforced with constraints
    StartI_ConstrEast = rhs.StartI_ConstrEast;
    EndI_ConstrEast = rhs.EndI_ConstrEast;
    StartJ_ConstrEast = rhs.StartJ_ConstrEast;
    EndJ_ConstrEast = rhs.EndJ_ConstrEast;

    // set range of cells which have constrained reconstruction if the North boundary condition is enforced with constraints
    StartI_ConstrNorth = rhs.StartI_ConstrNorth;
    EndI_ConstrNorth = rhs.EndI_ConstrNorth;
    StartJ_ConstrNorth = rhs.StartJ_ConstrNorth;
    EndJ_ConstrNorth = rhs.EndJ_ConstrNorth;

    // set range of cells which have constrained reconstruction if the South boundary condition is enforced with constraints
    StartI_ConstrSouth = rhs.StartI_ConstrSouth;
    EndI_ConstrSouth = rhs.EndI_ConstrSouth;
    StartJ_ConstrSouth = rhs.StartJ_ConstrSouth;
    EndJ_ConstrSouth = rhs.EndJ_ConstrSouth;

    // set range of cells which have SI calculation modified if the West boundary condition is enforced with constraints
    StartJ_SI_ConstrWest = rhs.StartJ_SI_ConstrWest;
    EndJ_SI_ConstrWest = rhs.EndJ_SI_ConstrWest;

    // set range of cells which have SI calculation modified if the East boundary condition is enforced with constraints
    StartJ_SI_ConstrEast = rhs.StartJ_SI_ConstrEast;
    EndJ_SI_ConstrEast = rhs.EndJ_SI_ConstrEast;

    // set range of cells which have SI calculation modified if the North boundary condition is enforced with constraints
    StartI_SI_ConstrNorth = rhs.StartI_SI_ConstrNorth;
    EndI_SI_ConstrNorth = rhs.EndI_SI_ConstrNorth;

    // set range of cells which have SI calculation modified if the South boundary condition is enforced with constraints
    StartI_SI_ConstrSouth = rhs.StartI_SI_ConstrSouth;
    EndI_SI_ConstrSouth = rhs.EndI_SI_ConstrSouth;

    // set the type of reconstruction required by the geometry setup
    _constrained_block_reconstruction = rhs._constrained_block_reconstruction;
    _constrained_WEST_reconstruction = rhs._constrained_WEST_reconstruction;
    _constrained_EAST_reconstruction = rhs._constrained_EAST_reconstruction;
    _constrained_NORTH_reconstruction = rhs._constrained_NORTH_reconstruction;
    _constrained_SOUTH_reconstruction = rhs._constrained_SOUTH_reconstruction;

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
    }

    delete [] TD; TD = NULL;
    delete [] SI; SI = NULL;
    delete [] LimitedCell; LimitedCell = NULL;
    delete [] PreviousLimitedCell; PreviousLimitedCell = NULL;

    if (_allocated_psinv){
      delete [] CENO_LHS; CENO_LHS = NULL;
      delete [] CENO_Geometric_Weights; CENO_Geometric_Weights = NULL;
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
    StartI_SI = 0; EndI_SI = 0;
    StartJ_SI = 0; EndJ_SI = 0;
    StartI_LPWL = 0; EndI_LPWL = 0;
    StartJ_LPWL = 0; EndJ_LPWL = 0;
    StartI_ConstrWest = 0; EndI_ConstrWest = 0;
    StartJ_ConstrWest = 0; EndJ_ConstrWest = 0;
    StartI_ConstrEast = 0; EndI_ConstrEast = 0;
    StartJ_ConstrEast = 0; EndJ_ConstrEast = 0;
    StartI_ConstrNorth = 0; EndI_ConstrNorth = 0;
    StartJ_ConstrNorth = 0; EndJ_ConstrNorth = 0;
    StartI_ConstrSouth = 0; EndI_ConstrSouth = 0;
    StartJ_ConstrSouth = 0; EndJ_ConstrSouth = 0;
    StartJ_SI_ConstrWest = 0;  EndJ_SI_ConstrWest = 0;
    StartJ_SI_ConstrEast = 0;  EndJ_SI_ConstrEast = 0;
    StartI_SI_ConstrNorth = 0; EndI_SI_ConstrNorth = 0;
    StartI_SI_ConstrSouth = 0; EndI_SI_ConstrSouth = 0;

    // Separate the high-order object from the associated geometry
    Geom = NULL;
    // Reset grid observers
    ObserverInteriorCellGeometryState = 0;
    ObserverGhostCellGeometryState = 0;
    ObserverCornerGhostCellGeometryState = 0;
    _constrained_block_reconstruction = false;
    _constrained_WEST_reconstruction  = false;
    _constrained_EAST_reconstruction  = false;
    _constrained_NORTH_reconstruction = false;
    _constrained_SOUTH_reconstruction = false;

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
    StartI_SI = 0; EndI_SI = 0;
    StartJ_SI = 0; EndJ_SI = 0;
    StartI_LPWL = 0; EndI_LPWL = 0;
    StartJ_LPWL = 0; EndJ_LPWL = 0;
    StartI_ConstrWest = 0; EndI_ConstrWest = 0;
    StartJ_ConstrWest = 0; EndJ_ConstrWest = 0;
    StartI_ConstrEast = 0; EndI_ConstrEast = 0;
    StartJ_ConstrEast = 0; EndJ_ConstrEast = 0;
    StartI_ConstrNorth = 0; EndI_ConstrNorth = 0;
    StartJ_ConstrNorth = 0; EndJ_ConstrNorth = 0;
    StartI_ConstrSouth = 0; EndI_ConstrSouth = 0;
    StartJ_ConstrSouth = 0; EndJ_ConstrSouth = 0;
    StartJ_SI_ConstrWest = 0;  EndJ_SI_ConstrWest = 0;
    StartJ_SI_ConstrEast = 0;  EndJ_SI_ConstrEast = 0;
    StartI_SI_ConstrNorth = 0; EndI_SI_ConstrNorth = 0;
    StartI_SI_ConstrSouth = 0; EndI_SI_ConstrSouth = 0;

    _allocated_cells = false;
    _allocated_psinv = false;
  }//endif
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

  // Determine the range of cells in which NO constrained reconstruction is performed
  SetRangeOfQuadCellsWithoutConstrainedReconstruction();

  // Determine the range of cells in which !constrained! reconstruction is performed
  SetRangeOfQuadCellsWithConstrainedReconstruction();

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

  // Compute the smoothness indicator for cells that use the central stencil
  for ( j  = StartJ_SI ; j <= EndJ_SI ; ++j ) {
    for ( i = StartI_SI ; i <= EndI_SI ; ++i ) {

      // Set the supporting stencil
      SetSmoothnessIndicatorStencil(i,j);

      // Evaluate the Smoothness Indicator for the current cell for all solution state variables
      ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
				 i, j, i_index, j_index, StencilSize_SmoothnessIndicator);

      // Check the smoothness condition
      AssessInterpolantsSmoothness(i,j);
 
    } /* endfor(i) */
  } /* endfor(j) */

  // Check whether smoothness indicator calculation with modified stencil is required anywhere in the block.
  if ( !_constrained_block_reconstruction ){
    // No need to perform any calculation with modified stencil
    return;
  }


  /***************************************************************************
   *    Perform smoothness indicator calculation with modified stencil       *
   **************************************************************************/

  // === Check WEST boundary
  if (_constrained_WEST_reconstruction){
    for ( j  = StartJ_SI_ConstrWest ; j <= EndJ_SI_ConstrWest ; ++j ) {
      for ( i = ICl ; i < StartI_SI ; ++i ) {

	// Set the supporting stencil
	SetSmoothnessIndicatorStencilForConstrainedReconstruction(i,j);

	// Evaluate the Smoothness Indicator for the current cell for all solution state variables
	ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
				   i, j, i_index_ave, j_index_ave, i_index_ave.size());

	// Check the smoothness condition
	AssessInterpolantsSmoothness(i,j);
 
      } /* endfor(i) */
    } /* endfor(j) */
  } // endfor (_constrained_WEST_reconstruction)

  // === Check EAST boundary
  if (_constrained_EAST_reconstruction){
    for ( j  = StartJ_SI_ConstrEast ; j <= EndJ_SI_ConstrEast ; ++j ) {
      for ( i = EndI_SI + 1 ; i <= ICu ; ++i ) {

	// Set the supporting stencil
	SetSmoothnessIndicatorStencilForConstrainedReconstruction(i,j);

	// Evaluate the Smoothness Indicator for the current cell for all solution state variables
	ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
				   i, j, i_index_ave, j_index_ave, i_index_ave.size());

	// Check the smoothness condition
	AssessInterpolantsSmoothness(i,j);
 
      } /* endfor(i) */
    } /* endfor(j) */
  } // endfor (_constrained_EAST_reconstruction)

  // === Check NORTH boundary
  if (_constrained_NORTH_reconstruction){
    for ( j  = EndJ_SI + 1 ; j <= JCu ; ++j ) {
      for ( i = StartI_SI_ConstrNorth ; i <= EndI_SI_ConstrNorth ; ++i ) {

	// Set the supporting stencil
	SetSmoothnessIndicatorStencilForConstrainedReconstruction(i,j);

	// Evaluate the Smoothness Indicator for the current cell for all solution state variables
	ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
				   i, j, i_index_ave, j_index_ave, i_index_ave.size());

	// Check the smoothness condition
	AssessInterpolantsSmoothness(i,j);
 
      } /* endfor(i) */
    } /* endfor(j) */
  } // endfor (_constrained_NORTH_reconstruction)

  // === Check SOUTH boundary
  if (_constrained_SOUTH_reconstruction){
    for ( j  = JCl ; j < StartJ_SI ; ++j ) {
      for ( i = StartI_SI_ConstrSouth ; i <= EndI_SI_ConstrSouth ; ++i ) {

	// Set the supporting stencil
	SetSmoothnessIndicatorStencilForConstrainedReconstruction(i,j);

	// Evaluate the Smoothness Indicator for the current cell for all solution state variables
	ComputeSmoothnessIndicator(SolnBlk, ReconstructedSoln,
				   i, j, i_index_ave, j_index_ave, i_index_ave.size());

	// Check the smoothness condition
	AssessInterpolantsSmoothness(i,j);
 
      } /* endfor(i) */
    } /* endfor(j) */
  } // endfor (_constrained_SOUTH_reconstruction)

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
 * Set the first and last indexes in 'I' and 'J' directions
 * which correspond to cells for which NO high-order constrained 
 * reconstruction is performed.
 * This type of reconstruction is typically required in the proximity
 * of curved boundary.
 * Flag also the boundaries for which constrained reconstruction
 * is required.
 *
 * \note The associated grid is used to determine this information!
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetRangeOfQuadCellsWithoutConstrainedReconstruction(void) {
  
  // Initialize indexes as if no constrained reconstruction would exist.
  StartI = ICl - Nghost_HO; EndI = ICu + Nghost_HO;
  StartJ = JCl - Nghost_HO; EndJ = JCu + Nghost_HO;

  /* Initialize indexes for calculation of smoothness indicator
     (i.e. cells used for flux calculation) as if no constrained reconstruction would exit. */
  StartI_SI = ICl - 1; EndI_SI = ICu + 1;
  StartJ_SI = JCl - 1; EndJ_SI = JCu + 1;

  /* Initialize indexes between which the interpolants are checked for non-smoothness
     (i.e. cells used for flux calculation) as if no constrained reconstruction would exit. */
  StartI_LPWL = ICl - 1; EndI_LPWL = ICu + 1;
  StartJ_LPWL = JCl - 1; EndJ_LPWL = JCu + 1;

  // Reset affected flags
  _constrained_block_reconstruction = false;
  _constrained_WEST_reconstruction  = false; 
  _constrained_EAST_reconstruction  = false; 
  _constrained_NORTH_reconstruction = false; 
  _constrained_SOUTH_reconstruction = false; 
  
  /* Check for existence of constrained reconstruction
     at boundaries and modify the affected indexes accordingly. */
  // Check West boundary
  if (Geom->IsWestBoundaryReconstructionConstrained()) {
    // Set StartI
    StartI = ICl + rings;
    // Flag the boundary accordingly
    _constrained_WEST_reconstruction = true;
    // Set StartI_SI
    StartI_SI += 1 + rings_SI;
    // Set StartI_LPWL
    StartI_LPWL += 1;
  } 
  // Check East boundary
  if (Geom->IsEastBoundaryReconstructionConstrained()){
    // Set EndI
    EndI = ICu - rings;
    // Flag the boundary accordingly
    _constrained_EAST_reconstruction = true;    
    // Set EndI_SI
    EndI_SI -= 1 + rings_SI;
    // Set EndI_LPWL
    EndI_LPWL -= 1;
  } 
  // Check South boundary
  if (Geom->IsSouthBoundaryReconstructionConstrained()){
    // Set StartJ
    StartJ = JCl + rings;
    // Flag the boundary accordingly
    _constrained_SOUTH_reconstruction = true;
    // Set StartJ_SI
    StartJ_SI += 1 + rings_SI;
    // Set StartJ_LPWL
    StartJ_LPWL += 1;
  } 
  // Check North boundary
  if (Geom->IsNorthBoundaryReconstructionConstrained()){
    // Set EndJ
    EndJ = JCu - rings;
    // Flag the boundary accordingly
    _constrained_NORTH_reconstruction = true;    
    // Set EndJ_SI
    EndJ_SI -= 1 + rings_SI;
    // Set EndJ_LPWL
    EndJ_LPWL -= 1;
  }

  // Flag the block accordingly
  if (_constrained_WEST_reconstruction  || _constrained_EAST_reconstruction || 
      _constrained_SOUTH_reconstruction || _constrained_NORTH_reconstruction ){
    _constrained_block_reconstruction = true;
  }
}

/*! 
 * Set the indexes in 'I' and 'J' directions which keep 
 * track of how the constrained reconstructions should be performed.
 *
 * \note This routine should be called after SetRangeOfQuadCellsWithoutConstrainedReconstruction()
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetRangeOfQuadCellsWithConstrainedReconstruction(void) {
  
  /* Check for existence of constrained reconstruction
     at boundaries and modify the affected indexes accordingly. */
  // Check West boundary
  if (_constrained_WEST_reconstruction){
    
    StartI_ConstrWest = ICl;
    EndI_ConstrWest   = StartI - 1;

    StartJ_ConstrWest = StartJ;
    EndJ_ConstrWest   = EndJ;

    // SI calculation
    StartJ_SI_ConstrWest = StartJ_SI;
    EndJ_SI_ConstrWest   = EndJ_SI;
  } else {

    StartI_ConstrWest = 0;
    EndI_ConstrWest   = 0;
    StartJ_ConstrWest = 0;
    EndJ_ConstrWest   = 0;

    // SI calculation
    StartJ_SI_ConstrWest = 0;
    EndJ_SI_ConstrWest   = 0;    
  }

  // Check East boundary
  if (_constrained_EAST_reconstruction){
    
    StartI_ConstrEast = EndI+1;
    EndI_ConstrEast   = ICu;

    StartJ_ConstrEast = StartJ;
    EndJ_ConstrEast   = EndJ;

    // SI calculation
    StartJ_SI_ConstrEast = StartJ_SI;
    EndJ_SI_ConstrEast   = EndJ_SI;
  } else {

    StartI_ConstrEast = 0;
    EndI_ConstrEast   = 0;
    StartJ_ConstrEast = 0;
    EndJ_ConstrEast   = 0;

    // SI calculation
    StartJ_SI_ConstrEast = 0;
    EndJ_SI_ConstrEast   = 0;
  }

  // Check South boundary
  if (_constrained_SOUTH_reconstruction){

    StartJ_ConstrSouth = JCl;
    EndJ_ConstrSouth   = StartJ - 1;    

    // check West
    if (_constrained_WEST_reconstruction){
      StartI_ConstrSouth = ICl;

      // SI calculation
      StartI_SI_ConstrSouth = ICl;
    } else {
      StartI_ConstrSouth = StartI;

      // SI calculation
      StartI_SI_ConstrSouth = StartI_SI;
    }

    // check East
    if (_constrained_EAST_reconstruction){
      EndI_ConstrSouth   = ICu;

      // SI calculation
      EndI_SI_ConstrSouth = ICu;  
    } else {
      EndI_ConstrSouth   = EndI;

      // SI calculation
      EndI_SI_ConstrSouth = EndI_SI;
    }

  } else {

    StartI_ConstrSouth = 0;
    EndI_ConstrSouth   = 0;
    StartJ_ConstrSouth = 0;
    EndJ_ConstrSouth   = 0;

    // SI calculation
    StartI_SI_ConstrSouth = 0;
    EndI_SI_ConstrSouth   = 0;
  }

  // Check North boundary
  if (_constrained_NORTH_reconstruction){

    StartJ_ConstrNorth = EndJ + 1;
    EndJ_ConstrNorth   = JCu;

    // check West
    if (_constrained_WEST_reconstruction){
      StartI_ConstrNorth = ICl;

      // SI calculation
      StartI_SI_ConstrNorth = ICl;
    } else {
      StartI_ConstrNorth = StartI;

      // SI calculation
      StartI_SI_ConstrNorth = StartI_SI;
    }

    // check East
    if (_constrained_EAST_reconstruction){
      EndI_ConstrNorth   = ICu;

      // SI calculation
      EndI_SI_ConstrNorth = ICu;
    } else {
      EndI_ConstrNorth   = EndI;

      // SI calculation
      EndI_SI_ConstrNorth = EndI_SI;
    }

  } else {

    StartI_ConstrNorth = 0;
    EndI_ConstrNorth   = 0;
    StartJ_ConstrNorth = 0;
    EndJ_ConstrNorth   = 0;

    // SI calculation
    StartI_SI_ConstrNorth = 0;
    EndI_SI_ConstrNorth   = 0;
  }

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
 * the constrained reconstruction of cell (iCell,jCell).
 * Use the number of rings and the class variables caring information
 * about constrained boundaries to determine how far the stencil extends.
 * This routine DOES'T generate a central stencil!
 * The stencil is biased to the mesh interior but it doesn't extend
 * further than a central stencil.
 *
 * \param [out] i_index The i-index of the cells.
 * \param [out] j_index The j-index of the cells.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetConstrainedReconstructionStencil(const int &iCell, const int &jCell,
								  IndexType & i_index, IndexType & j_index) const{

  // Reset indexes
  i_index.clear();
  j_index.clear();
  
  int i,j, Imin, Imax, Jmin, Jmax;

  /* Set Imin, Imax, Jmin, Jmax for a central stencil. */
  Imin = iCell-rings; Imax = iCell+rings;
  Jmin = jCell-rings; Jmax = jCell+rings;

  if (CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_ADDITIONAL_APPROXIMATE_CONSTRAINTS == ON &&
      CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL == OFF) {

    // Additional equations come from approximate constraints.
    
    // Check WEST boundary
    if (_constrained_WEST_reconstruction && Imin < ICl){
      Imin = ICl;		/* limit Imin */
    } 
    
    // Check EAST boundary
    if (_constrained_EAST_reconstruction && Imax > ICu){
      Imax = ICu;		/* limit Imax */
    } 
    
    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction && Jmax > JCu){
      Jmax = JCu;		/* limit Jmax */
    } 
    
    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction && Jmin < JCl){
      Jmin = JCl;		/* limit Jmin */
    }

  } else {

    // Additional equations come from extending the stencil in the opposite direction of the constrained boundary

    // Check WEST boundary
    if (_constrained_WEST_reconstruction && Imin < ICl){
      Imin = ICl;		/* limit Imin */
      
      Imax += 1;		/* extend Imax */
    } 
    
    // Check EAST boundary
    if (_constrained_EAST_reconstruction && Imax > ICu){
      Imax = ICu;		/* limit Imax */

      Imin -= 1;                /* extend Imin */
    } 
    
    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction && Jmax > JCu){
      Jmax = JCu;		/* limit Jmax */

      Jmin -= 1;                /* extend Jmin */
    } 
    
    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction && Jmin < JCl){
      Jmin = JCl;		/* limit Jmin */

      Jmax += 1;                /* extend Jmax */
    }

  } // endif


  /* Form stencil */
  i_index.push_back(iCell);
  j_index.push_back(jCell);

  for (i=Imin; i<=Imax; ++i){
    for (j=Jmin; j<=Jmax; ++j){

      if (  i!=iCell || j!=jCell ){
	i_index.push_back(i);
	j_index.push_back(j);
      }//endif

    }// endfor
  }// endfor

}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the smoothness indicator calculation stencil of cell (iCell,jCell).
 * Use the number of rings smoothness indicator set in the class to 
 * determine how far the stencil extends.
 * This routine doesn't modify the stencil due to existence 
 * of curved boundaries.
 * The indexes are written in the i_index and j_index containers.
 * Because these containers might be larger than needed, 
 * the StencilSize_SmoothnessIndicator variable is used for storing 
 * the smoothness indicator stencil size.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetSmoothnessIndicatorStencil(const int &iCell, const int &jCell){

  // Call set central stencil
  SetCentralStencil(iCell,jCell,i_index,j_index,rings_SI,StencilSize_SmoothnessIndicator);
}

/*! 
 * Write the 'i' and 'j' indexes of the cells that are part of
 * the smoothness indicator calculation stencil of cell (iCell,jCell),
 * which is affected by the presence of a constrained boundary condition.
 * Use the number of rings smoothness indicator set in the class to 
 * determine how far the stencil extends.
 * The indexes are written in the i_index_ave and j_index_ave containers.
 *
 * \note The first position (i_index[0],j_index[0]) corresponds to (iCell,jCell).
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetSmoothnessIndicatorStencilForConstrainedReconstruction(const int &iCell,
											const int &jCell){

  // Reset indexes
  i_index_ave.clear();
  j_index_ave.clear();
  
  int i,j, Imin, Imax, Jmin, Jmax;

  /* Set Imin, Imax, Jmin, Jmax for a central stencil with the number of rings dictated by rings_SI. */
  Imin = iCell-rings_SI; Imax = iCell+rings_SI;
  Jmin = jCell-rings_SI; Jmax = jCell+rings_SI;

  if (CENO_Execution_Mode::CENO_CONSTRAINED_RECONSTRUCTION_WITH_EXTENDED_BIASED_STENCIL == OFF || 
      CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS ) {

    // Don't enlarge the stencil, only cut off the ghost cells associated with constrained boundaries.
    
    // Check WEST boundary
    if (_constrained_WEST_reconstruction && Imin < ICl){
      Imin = ICl;		/* limit Imin */
    } 
    
    // Check EAST boundary
    if (_constrained_EAST_reconstruction && Imax > ICu){
      Imax = ICu;		/* limit Imax */
    } 
    
    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction && Jmax > JCu){
      Jmax = JCu;		/* limit Jmax */
    } 
    
    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction && Jmin < JCl){
      Jmin = JCl;		/* limit Jmin */
    }

  } else {

    // Enlarge the stencil in the opposite direction of the constrained boundary,
    // to be consistent with the one used for the constrained reconstruction

    // Check WEST boundary
    if (_constrained_WEST_reconstruction && Imin < ICl){
      Imin = ICl;		/* limit Imin */
      
      Imax += 1;		/* extend Imax */
    } 
    
    // Check EAST boundary
    if (_constrained_EAST_reconstruction && Imax > ICu){
      Imax = ICu;		/* limit Imax */

      Imin -= 1;                /* extend Imin */
    } 
    
    // Check NORTH boundary
    if (_constrained_NORTH_reconstruction && Jmax > JCu){
      Jmax = JCu;		/* limit Jmax */

      Jmin -= 1;                /* extend Jmin */
    } 
    
    // Check SOUTH boundary
    if (_constrained_SOUTH_reconstruction && Jmin < JCl){
      Jmin = JCl;		/* limit Jmin */

      Jmax += 1;                /* extend Jmax */
    }

  } // endif


  /* Form stencil */
  i_index_ave.push_back(iCell);
  j_index_ave.push_back(jCell);

  for (i=Imin; i<=Imax; ++i){
    for (j=Jmin; j<=Jmax; ++j){

      if (  i!=iCell || j!=jCell ){
	i_index_ave.push_back(i);
	j_index_ave.push_back(j);
      }//endif

    }// endfor
  }// endfor

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
  Communicator.Bcast(buffer, buffer_size, MPI::DOUBLE, Source_CPU);

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
