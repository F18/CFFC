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
  typedef std::vector<double> DoubleArrayType;

  typedef double (ClassType::*MemberFunction_OneArgument_Type)(const double &);
  typedef double (ClassType::*MemberFunction_TwoArguments_Type)(const double &, const double &);
  typedef double (ClassType::*MemberFunction_TwoArguments_OneParameter_Type)(const double &, const unsigned );
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
  const short int & CellInadequateFitValue(const int & ii, const int & jj,
					   const int VarPosition) const { return LimitedCell[ii][jj][VarPosition-1];}
  //! Get the monotonicity flag for reconstruction of cell (ii,jj) and the variable stored in the position 'VarPosition'.
  short int & CellInadequateFitValue(const int & ii, const int & jj,
				     const int VarPosition){ return LimitedCell[ii][jj][VarPosition-1];}
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

  //! @name Cell geometry
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
  const short int & NghostHO(void) const { return Nghost_HO; }
  void ResetMonotonicityFlag(void);
  void InitializeMonotonicityVariables(const int & ii, const int & jj);
  void InitializeVariable(int ReconstructionOrder, GeometryType & Block,
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
  //! Evaluate the interpolant at a given location (X_Coord,Y_Coord) for all solution variables,
  //  using the reconstruction of cell (ii,jj)
  Soln_State SolutionAtCoordinates(const int & ii, const int & jj,
				   const double & X_Coord, const double & Y_Coord) const {
    return TD[ii][jj].ComputeSolutionFor(X_Coord - XCellCenter(ii,jj), Y_Coord - YCellCenter(ii,jj));
  }
  //@}

  /*! @brief Integrate over the domain of the geometry associated with this high-order solution  */
  template<typename FO, class ReturnType>
  ReturnType IntegrateOverTheCell(const int &ii, const int &jj, const FO FuncObj,
				  const int & digits, ReturnType _dummy_param) const;
  
  //! @name Reconstructions:
  //@{
  /*! @brief Compute the unlimited high-order solution reconstruction.  */
  template<class Soln_Block_Type>
  void ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk);

  /*! @brief Compute the pseudo-inverse corresponding to the unlimited high-order solution reconstruction.  */
  template<class Soln_Block_Type>
  void ComputeReconstructionPseudoInverse(void);

  /*! @brief Compute the second (low-order) reconstruction in the CENO algorithm.  */
  template<class Soln_Block_Type>
  void ComputeLowOrderReconstruction(Soln_Block_Type &SolnBlk,
				     const int &Limiter);
  //@}

  //! @name CENO Analysis:
  //@{
  template<class Soln_Block_Type>
  void ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk);
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

  /* Friend functions */
  friend ostream & operator<< <Soln_State> (ostream & os, const HighOrder2D<Soln_State> & Obj);
  friend istream & operator>> <Soln_State> (istream & os, HighOrder2D<Soln_State> & Obj);

protected:
  
private:
  int Ni;		       //!< Number of high-order objects in i-direction
  int Nj;		       //!< Number of high-order objects in j-direction
  int Ng; 		       //!< Number of block ghost cells 
  int ICl,		       //!< Index of first interior cell in i-direction 
    ICu,		       //!< Index of last interior cell in i-direction 
    JCl,		       //!< Index of first interior cell in j-direction 
    JCu;                       //!< Index of last interior cell in j-direction 
  int OrderOfReconstruction;   //!< The order of reconstruction of the high-order object.
  short int Nghost_HO;	       //!< Number of ghost cells in which high-order reconstruction is performed. 

  bool _allocated_block;       //!< Flag indicating if the containers at block level has been allocated or not.
  bool _allocated_cells;       //!< Flag indicating if the containers at cell level has been allocated or not.
  bool _allocated_psinv;       //!< Flag indicating if the pseudo-inverse related containers have been allocated or not. 
  short int _si_calculation;   /*!< Flag to store the method for smoothness indicator calculation. Set to the 
				* same value as CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS. */

  DerivativesContainer **TD;   //!< High-order TaylorDerivatives
  DoubleArrayType **SI;        //!< The values of the smoothness indicator calculated for each reconstructed variable
  FlagType **LimitedCell;      //!< Monotonicity flag: Values --> OFF - high-order reconstruction,
                               //                                  ON - limited linear reconstruction
  int rings;                   //!< Number of rings used to generate the reconstruction stencil
  bool _calculated_psinv;      //!< Flag to indicate whether the pseudo-inverse has been already computed or not.

  DenseMatrix **CENO_LHS;      //!< Storage for the pseudo-inverse of the LHS term in the CENO reconstruction.
  DoubleArrayType **CENO_Geometric_Weights;   //!< Storage for the geometric weights used in CENO reconstruction.

  GeometryType* Geom;          //!< Pointer to the associated geometry which is a 2D mesh block

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
  //! Get the number of variables in the solution state
  int NumberOfVariables(void) const {return Soln_State::NumVar(); }

  //! Allocate memory at the cell level based on the order of reconstruction
  void allocate_CellMemory(const int &ReconstructionOrder, const bool &_pseudo_inverse_allocation_);

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
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  TD(NULL), SI(NULL), LimitedCell(NULL),
  rings(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL),
  Geom(NULL), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS)
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
  _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
  TD(NULL), SI(NULL), LimitedCell(NULL), 
  rings(0), _calculated_psinv(false),
  CENO_LHS(NULL), CENO_Geometric_Weights(NULL), 
  Geom(&Block), _si_calculation(CENO_Execution_Mode::CENO_SMOOTHNESS_INDICATOR_COMPUTATION_WITH_ONLY_FIRST_NEIGHBOURS)
{
  // Use the grid to get the number of interior block cells and the number of available ghost cells
  allocate(ICu_Grid()-ICl_Grid()+1,
	   JCu_Grid()-JCl_Grid()+1,
	   Nghost_Grid(),
	   _pseudo_inverse_allocation_,
	   ReconstructionOrder);
}

//! Copy constructor 
template<class SOLN_STATE> inline
HighOrder2D<SOLN_STATE>::HighOrder2D(const HighOrder2D<SOLN_STATE> & rhs)
  : Ni(0), Nj(0), Ng(0),
    ICl(0), ICu(0), JCl(0), JCu(0),
    OrderOfReconstruction(-1), Nghost_HO(0),
    _allocated_block(false), _allocated_cells(false), _allocated_psinv(false),
    TD(NULL), SI(NULL), LimitedCell(NULL),
    rings(0), _calculated_psinv(false),
    CENO_LHS(NULL), CENO_Geometric_Weights(NULL),
    Geom(rhs.Geom), _si_calculation(rhs._si_calculation)
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
	  LimitedCell[i][j] = LimitedCell[i][j];

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
    
    // check if the rhs has cell memory allocated
    if (rhs._allocated_cells){

      // copy the cell containers
      for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
	for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i) {
	  TD[i][j] = rhs.TD[i][j];
	  SI[i][j] = rhs.SI[i][j];
	  LimitedCell[i][j] = LimitedCell[i][j];

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
      if (_pseudo_inverse_allocation_){
	CENO_LHS = new DenseMatrix* [Ni];
	CENO_Geometric_Weights = new DoubleArrayType* [Ni];
      }

      for (i = 0; i < Ni ; ++i){
	TD[i] = new DerivativesContainer [Nj];
	SI[i] = new DoubleArrayType [Nj];
	LimitedCell[i] = new FlagType [Nj];
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

  // Set the new reconstruction order
  OrderOfReconstruction = ReconstructionOrder;

  // Set the Nghost_HO based on the OrderOfReconstruction
  Nghost_HO = getNghostHighOrder();

  // Set the number of neighbour rings based on the OrderOfReconstruction
  SetRings();

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
	CENO_LHS[i][j].newsize(getStencilSize() - 1, TD[i][j].size() - 1);
	CENO_Geometric_Weights[i][j].assign(getStencilSize(), 0.0);
      }

    } /* endfor */
  }/* endfor */

  // Confirm allocation
  _allocated_cells = true;
  if (_pseudo_inverse_allocation_){
    _allocated_psinv = true;
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

      if (_allocated_psinv){
	delete [] CENO_LHS[i]; CENO_LHS[i] = NULL; // deallocate pseudo-inverse
	delete [] CENO_Geometric_Weights[i]; CENO_Geometric_Weights[i] = NULL; // deallocate geometric weights
      }
    }

    delete [] TD; TD = NULL;
    delete [] SI; SI = NULL;
    delete [] LimitedCell; LimitedCell = NULL;

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

    // Separate the high-order object from the associated geometry
    Geom = NULL;

    // Confirm the deallocation
    _allocated_block = false;
    _allocated_cells = false;
    _allocated_psinv = false;
    _calculated_psinv = false;
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
	
	// deallocate SmoothnessIndicator
	SI[ii][jj].clear();

	// deallocate LHS matrix and GeometricWeights
	if (_allocated_psinv){
	  CENO_LHS[ii][jj].newsize(0,0);
	  CENO_Geometric_Weights[ii][jj].clear();
	}

      } /* endfor */
    } /* endfor */

    // reset flags 
    _calculated_psinv = false;

    // Confirm the deallocation
    OrderOfReconstruction = -1;
    Nghost_HO = 0;
    rings = 0;

    _allocated_cells = false;
    _allocated_psinv = false;
  }//endif
}

//! Initialize the high-order object
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::InitializeVariable(int ReconstructionOrder, GeometryType & Block,
						 const bool &_pseudo_inverse_allocation_){

  // Allocate (re-allocate) memory for the high-order object.
  allocate(Block.ICu-Block.ICl+1,
	   Block.JCu-Block.JCl+1,
	   Block.Nghost,
	   _pseudo_inverse_allocation_,
	   ReconstructionOrder);

  // Associate geometry
  SetGeometryPointer(Block);

  // Compute the pseudo-inverse if required
  if (_pseudo_inverse_allocation_){
    // add the pseudo-inverse calculation here
  }
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
 */
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::SetRings(void){
  rings = getNumberOfRings(getTaylorDerivativesSize());
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

//! Reset the monotonicity flag throughout the block.
template<class SOLN_STATE> inline
void HighOrder2D<SOLN_STATE>::ResetMonotonicityFlag(void){

  int i,j,k;

  for (j  = JCl-Nghost_HO ; j <= JCu+Nghost_HO ; ++j ) {
    for ( i = ICl-Nghost_HO ; i <= ICu+Nghost_HO ; ++i ) {    
      for ( k = 0; k <= NumberOfVariables() - 1; ++k){
	LimitedCell[i][j][k] = OFF; /* reset the flags to OFF (smooth solution)*/
      } /* endfor */
    } /* endfor */
  } /* endfor */

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
ComputeUnlimitedSolutionReconstruction(Soln_Block_Type &SolnBlk){

#if 0
  vector<int> i_index(getStencilSize()); 
  string msg;

  switch(ReconstructionMethod){
  case RECONSTRUCTION_CENO:
    // Make Stencil
    MakeReconstructionStencil(Rings(),iCell,i_index);

    // Solve reconstruction for the current cell
    kExact_Reconstruction(*this,SolnBlk,AccessToHighOrderVar,iCell,i_index,getStencilSize(),NumberOfTaylorDerivatives());

    // Reset the CellInadequateFit flag & the limiter value in the Taylor derivatives container
    for (int i = 0; i <= NumberOfVariables() - 1; ++i){
      LimitedCell[i] = OFF; /* reset the flags to OFF (smooth solution)*/
    }
    TaylorDeriv().ResetLimiter();
    
    break;
    
  default:
    throw runtime_error("HighOrder2D ERROR: Unknown specified reconstruction method");

  }
#endif
}

// ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type> inline
void HighOrder2D<SOLN_STATE>::ComputeReconstructionPseudoInverse(void){

#if 0

  if (CENO_Execution_Mode::USE_CENO_ALGORITHM && 
      CENO_Execution_Mode::CENO_SPEED_EFFICIENT && 
      _calculated_psinv == OFF){

    // == Check if the reconstruction polynomial is piecewise constant
    if (NumberOfTaylorDerivatives() == 1){
      // Set properly the _calculated_psinv
      _calculated_psinv = ON;
      return;
    }

    vector<int> i_index(getStencilSize()); 

    // Make Stencil
    MakeReconstructionStencil(Rings(),iCell,i_index);
    
    // Form the left-hand-side (LHS) term for the current cell
    kExact_Reconstruction_Compute_LHS(*this,SolnBlk,AccessToHighOrderVar,
				      iCell,i_index,getStencilSize(),NumberOfTaylorDerivatives());
    
    // Compute the pseudo-inverse and override the LHS term
    LHS().pseudo_inverse_override();

    // Reset the CellInadequateFit flag & the limiter value in the Taylor derivatives container
    for (int i = 0; i <= NumberOfVariables() - 1; ++i){
      LimitedCell[i] = OFF; /* reset the flags to OFF (smooth solution)*/
    }
    TaylorDeriv().ResetLimiter();

    // Set properly the _calculated_psinv
    _calculated_psinv = ON;
  }
#endif

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
void HighOrder2D<SOLN_STATE>::ComputeLowOrderReconstruction(Soln_Block_Type &SolnBlk,
							    const int &Limiter){

#if 0
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
	
  for ( n = 1 ; n <= NumberOfVariables() ; ++n ) {

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
#endif
}

// ComputeUnlimitedSolutionReconstruction()
/*! 
 * Compute the unlimited high-order reconstruction for 
 * the computational cell iCell, using information provided by
 * the SolnBlk domain and the 'ReconstructionMethod' algorithm.
 */
template<class SOLN_STATE>
template<class Soln_Block_Type>
void HighOrder2D<SOLN_STATE>::ComputeSmoothnessIndicator(Soln_Block_Type &SolnBlk){

#if 0

  static double SS_Regression, SS_Residual; // regression sum of squares, residual sum of squares
  static double MeanSolution, alpha;
  static int DOF; 		 /* degrees of freedom */
  static double AdjustmentCoeff; /* adjustment coefficient */
  static double Temp, DeltaTol;	
  static int parameter, cell, ComputeSI;
  int _StencilSize_(getStencilSize());
  vector<int> i_index(_StencilSize_); 

  // Make Stencil
  MakeReconstructionStencil(Rings(),iCell,i_index);

  /* Compute the CENO smoothness indicator for the current cell */

  /* Initialize the static variables */
  DOF = NumberOfTaylorDerivatives();
  AdjustmentCoeff = (_StencilSize_ - DOF)/(DOF - 1.0);

  for (parameter=1; parameter<=NumberOfVariables(); ++parameter){

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
	  throw runtime_error("ComputeSmoothnessIndicator() 2D ERROR: CENO_CONSIDER_WEIGHTS only works with CENO_SPEED_EFFICIENT");
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
#endif

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


#if 0
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

#endif

/*! 
 * Output the current object to the 
 * provided output stream.
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Output_Object(ostream & out_file) const {

  int i,j;

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
  short int _si_calc_;
  int _Ni_, _Nj_, _Ng_, ReconstructionOrder;
  int i,j;

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
 * \param Block_Geometry the grid block to which the
 *                   high-order object is associated.
 *
 * \todo Implement this subroutine!
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Broadcast_HighOrder_Data(GeometryType & Block_Geometry){

#ifdef _MPI_VERSION

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
 *
 * \todo Implement this subroutine!
 */
template<class SOLN_STATE>
void HighOrder2D<SOLN_STATE>::Broadcast_HighOrder_Data(MPI::Intracomm &Communicator, 
						       const int &Source_CPU,
						       GeometryType & Block_Geometry){
  
  
}
#endif

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

#if 0

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
    for (i = ICl - ((SolnBlk[0].*AccessToHighOrderVar)().Rings() + 1);
	 i<= ICu + ((SolnBlk[0].*AccessToHighOrderVar)().Rings() + 1);
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


#endif // _HIGHORDER_2D_INCLUDED
