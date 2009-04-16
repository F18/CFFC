/*!\file HighOrder_MultiBlock.h
  \brief Header file implementing HighOrder_Multi_Block class. */

#ifndef _HIGHORDER_MULTIBLOCK_INCLUDED
#define _HIGHORDER_MULTIBLOCK_INCLUDED

/* Include required C++ libraries. */
// None

/* Using std namespace functions */
// None

/* Include CFFC header files */
#ifndef _HIGHORDER_INCLUDED
#include "HighOrder.h"	      // Include high-order class header file
#endif //_HIGHORDER_INCLUDED

#ifndef _HIGHORDER_INPUT_INCLUDED
#include "HighOrderInput.h"   // Include high-order input header file
#endif //_HIGHORDER_INPUT_INCLUDED

/*!
 * \class HighOrder_Multi_Block
 *
 * \brief Collection of templated functions that operate
 *        on high-order variables of 3D arrays of 
 *        quadrilateral solution blocks.
 * \nosubgrouping
 *********************************************************/
class HighOrder_Multi_Block {
public:

  //! Create initial high-order variables
  template<typename HEXA_BLOCK>
  static int Create_Initial_HighOrder_Variables(HEXA_BLOCK *Solution_Block,
						AdaptiveBlock3D_List &LocalSolnBlockList,
						const int ReconstructionOrder){

    int i;

    //  if (HighOrder_Input::NumberOfHighOrderReconstructions != 0){

    for (i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
      if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
	// allocate high-order objects
	Solution_Block[i].allocate_HighOrder(ReconstructionOrder);
	// allocate high-order boundary conditions
	// Solution_Block[i].allocate_HighOrder_BoundaryConditions(); //RR: allocate_HighOrder_Boundary_Conditions()
     
      } /* endif */
   
    }  /* endfor */
 
    //  } //endif (NumberOfHighOrderReconstructions)

    return(0);
  }


  /*!
   * Perform high-order reconstruction for 
   * all used solution blocks in the provided list.
   * 
   * \param ReconstructedSoln member function of Quad_Soln_Block which returns the solution.
   */

  template <class HEXA_BLOCK>
  static void HighOrder_Reconstruction(HEXA_BLOCK *Solution_Block,
				       AdaptiveBlock3D_List &LocalSolnBlockList,
				       const int LimiterType){
    //
    //			       const class HEXA_BLOCK::Soln_State &
    //			       (HEXA_BLOCK::*ReconstructedSoln)(const int &,
    //								const int &,
    //								const int &) const){
    
    for (int i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
      if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
	// Reconstruct the solution in the current solution block.
	Solution_Block[i].HighOrderVariable.ComputeHighOrderSolutionReconstruction(Solution_Block[i],
										   LimiterType,
										   &HEXA_BLOCK::CellSolution);
      } /* endif */
    }  /* endfor */
  }

  /*!
   * Routine: Linear_Reconstruction
   *                                                      
   * Perform piecewise limited linear reconstruct for 
   * the solution of a 1D array of 3D hexahedral multi-block 
   * solution blocks.
   */
  template<class HEXA_BLOCK>
  static void Linear_Reconstruction(HEXA_BLOCK *Solution_Block,
				    AdaptiveBlock3D_List &LocalSolnBlockList,
				    const int LimiterType){
  
     for (int i = 0 ; i <= LocalSolnBlockList.Nblk-1 ; ++i ) {
      if (LocalSolnBlockList.Block[i].used == ADAPTIVEBLOCK3D_USED) {
	Solution_Block[i].Linear_Reconstruction_LeastSquares(LimiterType);
      } /* endif */
     }  /* endfor */
    
  }
  
  void Linear_Reconstruction(const int Limiter);
  
private:
  HighOrder_Multi_Block(void);   //!< Private default constructor
  HighOrder_Multi_Block(const HighOrder_Multi_Block&); //!< Private copy constructor
  HighOrder_Multi_Block& operator=(const HighOrder_Multi_Block&); //!< Private assignment operator

};

#endif
