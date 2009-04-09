/*!\file HighOrder.h
  \brief Header file defining the templated HighOrder class.*/

#ifndef _HIGHORDER_INCLUDED
#define _HIGHORDER_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
using std::ostream;
using std::istream;
using std::vector;

/* Include CFFC header files */
#include "../Math/Matrix.h"
#include "../Math/NumericalLibrary.h"
#include "TDContainer.h"
#include "CENO_ExecutionMode.h"
#include "CENO_Tolerances.h"
#include "../Grid/Grid3DHexaBlock.h"
//#include "ReconstructionHelpers.h"
//#include "Cauchy_BoundaryConditions.h" // --> RR: replace later
//#include "HighOrder2D_BlockBoundary.h" // --> RR: replace later

/*********************************
 * Declare the HighOrder   class *
 ********************************/
template <class SOLN_STATE> 
class HighOrder;

/************************************************
 *     Friend Functions : HighOrder             *
 ************************************************/
template<class SOLN_STATE>
ostream & operator<< (ostream & os, const HighOrder<SOLN_STATE> & Obj);

template<class SOLN_STATE>
istream & operator>> (istream & os, HighOrder<SOLN_STATE> & Obj);

//!  \class HighOrder
//   -----------------------------------------------------------------
/*!  
 *   \brief Templated class for high-order variables in 3D.
 *
 *///-----------------------------------------------------------------
template<class SOLN_STATE>
class HighOrder{

public:

  //! @name Defined public types:
  //@{
  typedef SOLN_STATE Soln_State;
  typedef HighOrder<Soln_State> ClassType;
  typedef Grid3D_Hexa_Block GeometryType;
  typedef TaylorDerivativesContainer<Soln_State> DerivativesContainer;
  typedef typename DerivativesContainer::Derivative  Derivative;
  typedef typename GeometryType::GeometricMoments GeometricMoments;
  typedef typename GeometryType::GeomMoment       GeomMoment;
  //! Type for monotonicity information of each reconstructed variable.
  typedef std::vector<int> FlagType;
  //! Type for smoothness indicator associated with each reconstructed variable.
  typedef std::vector<double> DoubleArrayType;
  //! Type for array of Vector3D
  typedef vector<Vector3D> Vector3DArray;
  //! Type for high-order Cauchy boundary conditions
  //  typedef Cauchy_BCs<Soln_State> BC_Type;  --> RR: include BC later
  //! Type for the array of high-order solution boundary conditions
  //  typedef vector<BC_Type *> BC_Type_Array; --> RR: include BC later
  //@}


  //! @name Constructors:
  //@{
//  HighOrder(void);		//!< Simple constructor 
//  HighOrder(int ReconstructionOrder, GeometryType & Block,
//	      const bool &_pseudo_inverse_allocation_ = false); //!< Advanced constructor
//  HighOrder( const HighOrder & rhs); //!< Copy constructor
  void allocate(const int &NC_IDir, const int &NC_JDir, const int &NC_KDir,
		const int &Nghost,
		const bool &_pseudo_inverse_allocation_,
		int ReconstructionOrder = -1);
  //@}
  

  //! @name Destructors:
  //@{
//  ~HighOrder(void){ deallocate(); }
//  void deallocate_CellMemory(void);
//  void deallocate(void);
//  void deallocate_ReconstructionTypeMap(void);
//  //@} 
  
//  HighOrder<Soln_State> & operator=(const HighOrder<Soln_State> & rhs); //!< Assignment operator



  //! @name Initialize container functions.
  //@{ 
  void InitializeVariable(int ReconstructionOrder, GeometryType & Block,
			  const bool &_pseudo_inverse_allocation_ = false);
  void InitializeBasicVariable(int ReconstructionOrder, GeometryType & Block,
			       const bool &_pseudo_inverse_allocation_ = false);
  //@}

protected:

private:

  //! @name Internal indexes:
  //@{
  int Ni;		       //!< Number of high-order objects in i-direction
  int Nj;		       //!< Number of high-order objects in j-direction
  int Nk;		       //!< Number of high-order objects in k-direction
  int Ng; 		       //!< Number of block ghost cells 
  int ICl,		       //!< Index of first interior cell in i-direction 
    ICu,		       //!< Index of last interior cell in i-direction 
    JCl,		       //!< Index of first interior cell in j-direction 
    JCu,                       //!< Index of last interior cell in j-direction 
    KCl,		       //!< Index of first interior cell in k-direction 
    KCu;                       //!< Index of last interior cell in k-direction 
  int OrderOfReconstruction;   //!< The order of reconstruction of the high-order object.
  int Nghost_HO;	       //!< Number of ghost cells in which high-order reconstruction is performed. 
  int StartI,		       //!< Index of the first cell in i-direction in which NO constrained reconstruction is performed. 
    EndI,                      //!< Index of the last cell in i-direction in which NO constrained reconstruction is performed. 
    StartJ,                    //!< Index of the first cell in j-direction in which NO constrained reconstruction is performed. 
    EndJ,                      //!< Index of the last cell in j-direction in which NO constrained reconstruction is performed. 
    StartK,                    //!< Index of the first cell in k-direction in which NO constrained reconstruction is performed. 
    EndK;                      //!< Index of the last cell in k-direction in which NO constrained reconstruction is performed. 
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


  //RR: Block Boundary stuff
//  HighOrder2D_BlockBoundary WestBnd,	  //!< Storage for information related to the West high-order block boundary
//    S_WestBnd,                  //!< Storage for information related to the South extension of West high-order block boundary
//    N_WestBnd,                  //!< Storage for information related to the North extension of West high-order block boundary
//    EastBnd,			//!< Storage for information related to the East high-order block boundary
//    S_EastBnd,			//!< Storage for information related to the South extension of East high-order block boundary
//    N_EastBnd,			//!< Storage for information related to the North extension of East high-order block boundary
//    NorthBnd,                   //!< Storage for information related to the North high-order block boundary
//    W_NorthBnd,                 //!< Storage for information related to the West extension of North high-order block boundary
//    E_NorthBnd,                 //!< Storage for information related to the East extension of North high-order block boundary
//    SouthBnd,                   //!< Storage for information related to the South high-order block boundary
//    W_SouthBnd,                 //!< Storage for information related to the West extension of South high-order block boundary
//    E_SouthBnd;                 //!< Storage for information related to the East extension of South high-order block boundary

};




/****************************************************
 * Implement the HighOrder2D class member functions *
 ***************************************************/


// allocate()
/*! Allocate memory for the high-order object.
 *
 * \param NC_Idir number of cells in i-direction
 * \param NC_Jdir number of cells in j-direction
 * \param NC_Jdir number of cells in k-direction
 * \param Nghost  number of block ghost cells
 * \param ReconstructionOrder the order of reconstruction for this high-order object.
 *        If this number is not specified the implicit value is -1, which corresponds to no memory allocation.
 * 
 * \todo Improve the memory allocation based on the runtime situations and reconstruction order.
 */
template<class SOLN_STATE> inline
void HighOrder<SOLN_STATE>::allocate(const int &NC_IDir,
				       const int &NC_JDir,
				       const int &NC_KDir,
				       const int &Nghost,
				       const bool &_pseudo_inverse_allocation_,
				       int ReconstructionOrder){
  // RR: empty function

}


/*!
 * Initialize the high-order object.
 */
template<class SOLN_STATE> inline
void HighOrder<SOLN_STATE>::InitializeVariable(int ReconstructionOrder, GeometryType & Block,
					       const bool &_pseudo_inverse_allocation_){

  // Initialize the basic functionality of the high-order object
  InitializeBasicVariable(ReconstructionOrder, Block, _pseudo_inverse_allocation_);

  // Compute the pseudo-inverse if required
  //  ComputeReconstructionPseudoInverse();
}

/*! 
 * Initialize the basic functionality of the high-order object 
 * (i.e. memory allocation, indexes etc.).
 */
template<class SOLN_STATE> inline
void HighOrder<SOLN_STATE>::InitializeBasicVariable(int ReconstructionOrder, GeometryType & Block,
						    const bool &_pseudo_inverse_allocation_){

  // Allocate (re-allocate) memory for the high-order object.
  allocate(Block.ICu-Block.ICl+1,
	   Block.JCu-Block.JCl+1,
	   Block.KCu-Block.KCl+1,
	   Block.Nghost,
	   _pseudo_inverse_allocation_,
	   ReconstructionOrder);

  // RR: Include these basic initializations later

//  // Associate geometry
//  SetGeometryPointer(Block);
//
//  // Set properties of block boundaries (i.e. which boundaries are constrained/unconstrained and which are opaque/transparent)
//  SetPropertiesHighOrderBlock();
//
//  // Build the reconstruction type map
//  BuildReconstructionTypeMap();
//
//  // Compute the smoothness indicator adjustment coefficient
//  /* It adjusts the value of the smoothness indicator based on the degree of the reconstruction */
//  AdjustmentCoeff = (getStencilSize() - NumberOfTaylorDerivatives() )/( NumberOfTaylorDerivatives() - 1.0);
}



/* ---------------------------------------------------------------------------------------------- 
 * =============== INCLUDE THE IMPLEMENTATION OF HIGH-ORDER RECONSTRUCTIONS ==================
 * ---------------------------------------------------------------------------------------------*/
//#include "HighOrder_Reconstructions.h" --> RR: include later

#endif // _HIGHORDER_INCLUDED
